/*                 D B _ D I S C O V E R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file librt/db_discovery.cpp
 *
 * Focused large-model hierarchy discovery.  Workers only read immutable
 * directory records and write disjoint result slots.  Reference counts and
 * callbacks are gathered on the calling thread after all workers have
 * finished, preserving the ordinary librt publication contract.
 */

#include "common.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <string>
#include <vector>

#include "bu/mapped_file.h"
#include "bu/parallel.h"
#include "raytrace.h"
#include "rt/db_discovery.h"
#include "rt/db5.h"
#include "rt/geom.h"

#include "librt_private.h"


struct rt_db_hierarchy_arc_native {
    struct directory *parent = RT_DIR_NULL;
    struct directory *child = RT_DIR_NULL;
    std::string child_name;
    db_op_t operation = DB_OP_UNION;
    bool matrix_valid = false;
    mat_t matrix = MAT_INIT_IDN;
};


struct rt_db_hierarchy {
    std::vector<rt_db_hierarchy_arc_native> arcs;
};


struct rt_db_discovery_task {
    struct directory *dp = RT_DIR_NULL;
    size_t serialized_bytes = 0;
};


struct rt_db_discovery_task_result {
    std::vector<rt_db_hierarchy_arc_native> arcs;
};

enum rt_db_discovery_comb_token {
    RT_DB_DISCOVERY_TOKEN_LEAF = 1,
    RT_DB_DISCOVERY_TOKEN_UNION = 2,
    RT_DB_DISCOVERY_TOKEN_INTERSECT = 3,
    RT_DB_DISCOVERY_TOKEN_SUBTRACT = 4,
    RT_DB_DISCOVERY_TOKEN_XOR = 5,
    RT_DB_DISCOVERY_TOKEN_NOT = 6
};


static size_t
rt_db_discovery_width_bytes(int width)
{
    if (width < DB5HDR_WIDTHCODE_8BIT ||
	width > DB5HDR_WIDTHCODE_64BIT)
	return 0;
    return static_cast<size_t>(1) << width;
}


/*
 * A v5 record gives us the next record's address, but not the address after
 * that.  Issuing one tiny read for every record makes this dependency chain
 * pathological on a cold device.  This reader retains a modest sequential
 * window: dense object tables are streamed in large transfers, while a large
 * primitive body is still skipped when the next record falls outside the
 * window.  Consequently I/O adapts to record density without retaining or
 * copying primitive bodies.
 */
class rt_db_discovery_file_reader {
public:
    rt_db_discovery_file_reader(FILE *input, b_off_t file_length) :
	file(input), eof(file_length), window(4ULL * 1024ULL * 1024ULL)
    {
    }

    bool read(b_off_t offset, void *buffer, size_t bytes)
    {
	if (!file || !buffer || !bytes || offset < 0 || offset > eof ||
		bytes > static_cast<size_t>(eof - offset))
	    return false;

	if (window_start >= 0 && offset >= window_start) {
	    const b_off_t displacement = offset - window_start;
	    if (displacement >= 0 &&
		static_cast<size_t>(displacement) <= window_bytes &&
		bytes <= window_bytes - static_cast<size_t>(displacement)) {
		std::memcpy(buffer,
		    window.data() + static_cast<size_t>(displacement), bytes);
		return true;
	    }
	}

	if (bytes <= window.size()) {
	    const size_t available = static_cast<size_t>(eof - offset);
	    const size_t fill = std::min(window.size(), available);
	    if (bu_fseek(file, offset, SEEK_SET) != 0 ||
		fread(window.data(), 1, fill, file) != fill)
		return false;
	    window_start = offset;
	    window_bytes = fill;
	    std::memcpy(buffer, window.data(), bytes);
	    return true;
	}

	if (bu_fseek(file, offset, SEEK_SET) != 0)
	    return false;
	return fread(buffer, 1, bytes, file) == bytes;
    }

private:
    FILE *file = nullptr;
    b_off_t eof = 0;
    b_off_t window_start = -1;
    size_t window_bytes = 0;
    std::vector<unsigned char> window;
};


/*
 * Build the ordinary v5 directory without copying primitive bodies.
 *
 * Each object is sampled at its prefix, which contains the record length,
 * name, and normally its attributes.  Oversized name/attribute fields use a
 * bounded exact read.  Advancing by the validated encoded object length lets
 * the next prefix validate the preceding decision without reading the body.
 * This is particularly important for editable databases: the legacy stdio
 * scan reads every byte because it must construct a complete raw record.
 */
int
db5_scan_directory_headers(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->dbi_filename)
	return -1;
    RT_CK_DBI(dbip);

    FILE *file = fopen(dbip->dbi_filename, "rb");
    if (!file)
	return -1;
    if (bu_fseek(file, 0, SEEK_END) != 0) {
	fclose(file);
	return -1;
    }
    const b_off_t eof = bu_ftell(file);
    rt_db_discovery_file_reader reader(file, eof);
    unsigned char database_header[8];
    if (eof < static_cast<b_off_t>(sizeof(database_header)) ||
	!reader.read(0, database_header, sizeof(database_header)) ||
	db5_header_is_valid(database_header) == 0) {
	fclose(file);
	return -1;
    }

    static const size_t prefix_capacity = 512;
    /* The byte preceding each record is the prior record's trailer.  Reading
     * it with the next prefix validates boundaries without another I/O per
     * object. */
    std::vector<unsigned char> prefix(prefix_capacity + 1);
    std::vector<unsigned char> name_storage;
    std::vector<unsigned char> attribute_storage;
    b_off_t offset = static_cast<b_off_t>(sizeof(database_header));
    size_t records = 0;
    while (offset < eof) {
	const size_t available = static_cast<size_t>(eof - offset);
	const size_t prefix_bytes = std::min(prefix_capacity, available);
	if (prefix_bytes < sizeof(struct db5_ondisk_header) + 1 ||
	    !reader.read(offset - 1, prefix.data(),
		prefix_bytes + 1) || prefix[0] != DB5HDR_MAGIC2) {
	    fclose(file);
	    return -1;
	}
	unsigned char *record_prefix = prefix.data() + 1;
	const struct db5_ondisk_header *header =
	    reinterpret_cast<const struct db5_ondisk_header *>(record_prefix);
	if (header->db5h_magic1 != DB5HDR_MAGIC1) {
	    fclose(file);
	    return -1;
	}

	struct db5_raw_internal raw;
	std::memset(&raw, 0, sizeof(raw));
	raw.magic = DB5_RAW_INTERNAL_MAGIC;
	raw.h_dli = header->db5h_hflags & DB5HDR_HFLAGS_DLI_MASK;
	raw.h_object_width =
	    (header->db5h_hflags & DB5HDR_HFLAGS_OBJECT_WIDTH_MASK) >>
	    DB5HDR_HFLAGS_OBJECT_WIDTH_SHIFT;
	raw.h_name_present =
	    header->db5h_hflags & DB5HDR_HFLAGS_NAME_PRESENT;
	raw.h_name_hidden =
	    header->db5h_hflags & DB5HDR_HFLAGS_HIDDEN_OBJECT;
	raw.h_name_width =
	    (header->db5h_hflags & DB5HDR_HFLAGS_NAME_WIDTH_MASK) >>
	    DB5HDR_HFLAGS_NAME_WIDTH_SHIFT;
	raw.a_present = header->db5h_aflags & DB5HDR_AFLAGS_PRESENT;
	raw.a_width =
	    (header->db5h_aflags & DB5HDR_AFLAGS_WIDTH_MASK) >>
	    DB5HDR_AFLAGS_WIDTH_SHIFT;
	raw.a_zzz = header->db5h_aflags & DB5HDR_AFLAGS_ZZZ_MASK;
	raw.b_present = header->db5h_bflags & DB5HDR_BFLAGS_PRESENT;
	raw.b_width =
	    (header->db5h_bflags & DB5HDR_BFLAGS_WIDTH_MASK) >>
	    DB5HDR_BFLAGS_WIDTH_SHIFT;
	raw.b_zzz = header->db5h_bflags & DB5HDR_BFLAGS_ZZZ_MASK;
	raw.major_type = header->db5h_major_type;
	raw.minor_type = header->db5h_minor_type;
	BU_EXTERNAL_INIT(&raw.name);
	BU_EXTERNAL_INIT(&raw.attributes);
	BU_EXTERNAL_INIT(&raw.body);

	const size_t object_width =
	    rt_db_discovery_width_bytes(raw.h_object_width);
	if (!object_width ||
	    sizeof(struct db5_ondisk_header) + object_width > prefix_bytes) {
	    fclose(file);
	    return -1;
	}
	size_t object_units = 0;
	(void)db5_decode_length(&object_units,
	    record_prefix + sizeof(struct db5_ondisk_header),
	    raw.h_object_width);
	if (!object_units || object_units > SIZE_MAX / 8) {
	    fclose(file);
	    return -1;
	}
	raw.object_length = object_units * 8;
	if (raw.object_length < sizeof(struct db5_ondisk_header) +
		object_width || raw.object_length > available) {
	    fclose(file);
	    return -1;
	}

	size_t cursor = sizeof(struct db5_ondisk_header) + object_width;
	auto read_length = [&](int width, size_t &value) -> bool {
	    const size_t bytes = rt_db_discovery_width_bytes(width);
	    if (!bytes || cursor > raw.object_length ||
		bytes > raw.object_length - cursor)
		return false;
	    unsigned char encoded[8];
	    const unsigned char *source = nullptr;
	    if (cursor + bytes <= prefix_bytes)
		source = record_prefix + cursor;
	    else {
		if (!reader.read(
			offset + static_cast<b_off_t>(cursor), encoded, bytes))
		    return false;
		source = encoded;
	    }
	    (void)db5_decode_length(&value, source, width);
	    cursor += bytes;
	    return true;
	};
	auto read_field = [&](size_t bytes,
		std::vector<unsigned char> &storage,
		struct bu_external &external) -> bool {
	    if (cursor > raw.object_length || bytes > raw.object_length - cursor)
		return false;
	    external.ext_nbytes = bytes;
	    if (!bytes) {
		external.ext_buf = nullptr;
		return true;
	    }
	    if (cursor + bytes <= prefix_bytes) {
		external.ext_buf = record_prefix + cursor;
	    } else {
		storage.resize(bytes);
		if (!reader.read(
			offset + static_cast<b_off_t>(cursor), storage.data(),
			bytes))
		    return false;
		external.ext_buf = storage.data();
	    }
	    cursor += bytes;
	    return true;
	};

	if (raw.h_name_present) {
	    size_t bytes = 0;
	    if (!read_length(raw.h_name_width, bytes) ||
		!read_field(bytes, name_storage, raw.name) || !bytes ||
		raw.name.ext_buf[bytes - 1] != '\0') {
		fclose(file);
		return -1;
	    }
	}
	if (raw.a_present) {
	    size_t bytes = 0;
	    if (!read_length(raw.a_width, bytes) ||
		!read_field(bytes, attribute_storage, raw.attributes)) {
		fclose(file);
		return -1;
	    }
	}

	if (raw.h_dli == DB5HDR_HFLAGS_DLI_FREE_STORAGE) {
	    rt_memfree(&dbip->i->dbi_freep, raw.object_length, offset);
	} else if (raw.h_dli != DB5HDR_HFLAGS_DLI_HEADER_OBJECT &&
	    raw.name.ext_buf) {
	    (void)db5_diradd(dbip, &raw, offset, nullptr);
	}
	offset += static_cast<b_off_t>(raw.object_length);
	records++;
    }

    unsigned char final_trailer = 0;
    const bool valid_end = offset == eof && eof > 0 &&
	reader.read(eof - 1, &final_trailer, 1) &&
	final_trailer == DB5HDR_MAGIC2;
    fclose(file);
    if (!valid_end)
	return -1;
    dbip->i->dbi_eof = eof;
    dbip->i->dbi_nrec = records;
    return 0;
}


static uint64_t
rt_db_discovery_elapsed_us(
    const std::chrono::steady_clock::time_point &begin)
{
    return static_cast<uint64_t>(
	std::chrono::duration_cast<std::chrono::microseconds>(
	    std::chrono::steady_clock::now() - begin).count());
}


class rt_db_discovery_mapped_view {
public:
    explicit rt_db_discovery_mapped_view(struct db_i *database) : dbip(database)
    {
	if (!dbip || !dbip->i || db_version(dbip) != 5 ||
	    dbip->i->dbi_mf || !dbip->dbi_filename)
	    return;

	file = bu_open_mapped_file(dbip->dbi_filename, "rt_db_discovery");
	if (!file || !file->buf || !file->buflen) {
	    if (file)
		bu_close_mapped_file(file);
	    file = nullptr;
	    return;
	}

	previous_inmem = dbip->i->dbi_inmem;
	dbip->i->dbi_mf = file;
	dbip->i->dbi_inmem = file->buf;
	installed = true;
    }

    ~rt_db_discovery_mapped_view()
    {
	if (!installed)
	    return;
	dbip->i->dbi_mf = nullptr;
	dbip->i->dbi_inmem = previous_inmem;
	bu_close_mapped_file(file);
	/* Discovery is a one-shot startup optimization.  Do not leave a multi-GB
	 * no-user mapping resident after the editable database resumes. */
	bu_free_mapped_files(0);
    }

    bool active() const { return installed; }

private:
    struct db_i *dbip = nullptr;
    struct bu_mapped_file *file = nullptr;
    void *previous_inmem = nullptr;
    bool installed = false;
};


static db_op_t
rt_db_discovery_public_op(int operation)
{
    if (operation == OP_SUBTRACT)
	return DB_OP_SUBTRACT;
    if (operation == OP_INTERSECT)
	return DB_OP_INTERSECT;
    return DB_OP_UNION;
}


static void
rt_db_discovery_collect_comb_tree(struct db_i *dbip,
	struct directory *parent, union tree *tree, int operation,
	std::vector<rt_db_hierarchy_arc_native> &arcs)
{
    if (!tree)
	return;

    RT_CK_TREE(tree);
    switch (tree->tr_op) {
	case OP_DB_LEAF: {
	    if (!tree->tr_l.tl_name)
		return;
	    rt_db_hierarchy_arc_native arc;
	    arc.parent = parent;
	    arc.child_name = tree->tr_l.tl_name;
	    arc.child = db_lookup(dbip, tree->tr_l.tl_name, LOOKUP_QUIET);
	    arc.operation = rt_db_discovery_public_op(operation);
	    if (tree->tr_l.tl_mat) {
		arc.matrix_valid = true;
		MAT_COPY(arc.matrix, tree->tr_l.tl_mat);
	    }
	    arcs.push_back(std::move(arc));
	    return;
	}
	case OP_UNION:
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_left, OP_UNION, arcs);
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_right, OP_UNION, arcs);
	    return;
	case OP_INTERSECT:
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_left, OP_UNION, arcs);
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_right, OP_INTERSECT, arcs);
	    return;
	case OP_SUBTRACT:
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_left, OP_UNION, arcs);
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_right, OP_SUBTRACT, arcs);
	    return;
	case OP_XOR:
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_left, OP_UNION, arcs);
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_right, operation, arcs);
	    return;
	case OP_NOT:
	case OP_GUARD:
	case OP_XNOP:
	    rt_db_discovery_collect_comb_tree(dbip, parent,
		tree->tr_b.tb_left, OP_UNION, arcs);
	    return;
	default:
	    bu_log("rt_db_discovery: skipping unsupported tree operator %d\n",
		tree->tr_op);
	    return;
    }
}


static bool
rt_db_discovery_external_view(const struct db_i *dbip,
	const struct directory *dp, struct bu_external &external)
{
    size_t bytes = 0;
    const unsigned char *record = db_external_view(dbip, dp, &bytes);
    if (!record || !bytes)
	return false;

    BU_EXTERNAL_INIT(&external);
    external.ext_buf = const_cast<unsigned char *>(record);
    external.ext_nbytes = bytes;
    return true;
}


static void
rt_db_discovery_decode_comb(struct db_i *dbip, struct directory *dp,
	std::vector<rt_db_hierarchy_arc_native> &arcs)
{
    if (db_version(dbip) != 5) {
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, dbip, nullptr) < 0)
	    return;
	struct rt_comb_internal *comb = intern.idb_type == ID_COMBINATION ?
	    static_cast<struct rt_comb_internal *>(intern.idb_ptr) : nullptr;
	if (comb)
	    rt_db_discovery_collect_comb_tree(dbip, dp, comb->tree,
		OP_UNION, arcs);
	rt_db_free_internal(&intern);
	return;
    }

    struct bu_external external;
    if (!rt_db_discovery_external_view(dbip, dp, external))
	return;

    struct db5_raw_internal raw;
    if (!db5_get_raw_internal_ptr(&raw, external.ext_buf) ||
	raw.major_type != DB5_MAJORTYPE_BRLCAD ||
	raw.minor_type != DB5_MINORTYPE_BRLCAD_COMBINATION ||
	!raw.body.ext_buf)
	return;

    /* Decode the compact v5 combination body directly.  Importing a 150k-leaf
     * root merely to walk and free its allocated union tree was both the
     * dominant serial task and a source of allocator contention.  Discovery
     * needs only names, matrices, and effective boolean operations. */
    const unsigned char *cursor = raw.body.ext_buf;
    const unsigned char *end = raw.body.ext_buf + raw.body.ext_nbytes;
    if (cursor >= end)
	return;
    const int width = *cursor++;
    const size_t width_bytes = rt_db_discovery_width_bytes(width);
    if (!width_bytes)
	return;
    auto decode_length = [&](size_t &value) -> bool {
	if (static_cast<size_t>(end - cursor) < width_bytes)
	    return false;
	(void)db5_decode_length(&value, cursor, width);
	cursor += width_bytes;
	return true;
    };
    size_t matrix_count = 0;
    size_t leaf_count = 0;
    size_t leaf_bytes = 0;
    size_t expression_bytes = 0;
    size_t max_stack_depth = 0;
    if (!decode_length(matrix_count) || !decode_length(leaf_count) ||
	!decode_length(leaf_bytes) || !decode_length(expression_bytes) ||
	!decode_length(max_stack_depth))
	return;
    if (matrix_count > SIZE_MAX /
	(ELEMENTS_PER_MAT * SIZEOF_NETWORK_DOUBLE))
	return;
    const size_t matrix_bytes = matrix_count *
	ELEMENTS_PER_MAT * SIZEOF_NETWORK_DOUBLE;
    if (static_cast<size_t>(end - cursor) < matrix_bytes)
	return;
    const unsigned char *matrices = cursor;
    cursor += matrix_bytes;
    if (leaf_bytes > static_cast<size_t>(end - cursor))
	return;
    const unsigned char *leaf_cursor = cursor;
    const unsigned char *leaf_end = cursor + leaf_bytes;
    cursor = leaf_end;
    if (expression_bytes > static_cast<size_t>(end - cursor))
	return;
    const unsigned char *expression = cursor;

    struct expression_node {
	int operation = OP_DB_LEAF;
	size_t left = SIZE_MAX;
	size_t right = SIZE_MAX;
	size_t arc = SIZE_MAX;
    };
    std::vector<expression_node> nodes;
    std::vector<size_t> stack;
    if (expression_bytes) {
	nodes.reserve(expression_bytes);
	stack.reserve(std::max<size_t>(1, max_stack_depth));
    }
    arcs.reserve(arcs.size() + leaf_count);

    auto decode_leaf = [&]() -> size_t {
	if (leaf_cursor >= leaf_end)
	    return SIZE_MAX;
	const void *terminator = std::memchr(leaf_cursor, '\0',
	    static_cast<size_t>(leaf_end - leaf_cursor));
	if (!terminator)
	    return SIZE_MAX;
	const char *name = reinterpret_cast<const char *>(leaf_cursor);
	leaf_cursor = static_cast<const unsigned char *>(terminator) + 1;
	if (static_cast<size_t>(leaf_end - leaf_cursor) < width_bytes)
	    return SIZE_MAX;
	size_t matrix_index = 0;
	(void)db5_decode_signed(&matrix_index, leaf_cursor, width);
	leaf_cursor += width_bytes;

	rt_db_hierarchy_arc_native arc;
	arc.parent = dp;
	arc.child_name = name;
	arc.child = db_lookup(dbip, name, LOOKUP_QUIET);
	if (static_cast<ssize_t>(matrix_index) >= 0) {
	    if (matrix_index >= matrix_count)
		return SIZE_MAX;
	    double decoded[ELEMENTS_PER_MAT];
	    bu_cv_ntohd(reinterpret_cast<unsigned char *>(decoded),
		matrices + matrix_index * ELEMENTS_PER_MAT *
		    SIZEOF_NETWORK_DOUBLE, ELEMENTS_PER_MAT);
	    MAT_COPY(arc.matrix, decoded);
	    arc.matrix_valid = true;
	}
	const size_t arc_index = arcs.size();
	arcs.push_back(std::move(arc));
	return arc_index;
    };

    if (!expression_bytes) {
	for (size_t i = 0; i < leaf_count; ++i) {
	    if (decode_leaf() == SIZE_MAX) {
		arcs.clear();
		return;
	    }
	}
	if (leaf_cursor != leaf_end)
	    arcs.clear();
	return;
    }

    for (size_t i = 0; i < expression_bytes; ++i) {
	expression_node node;
	node.operation = expression[i];
	if (node.operation == RT_DB_DISCOVERY_TOKEN_LEAF) {
	    node.arc = decode_leaf();
	    if (node.arc == SIZE_MAX) {
		arcs.clear();
		return;
	    }
	} else if (node.operation == RT_DB_DISCOVERY_TOKEN_UNION ||
	    node.operation == RT_DB_DISCOVERY_TOKEN_INTERSECT ||
	    node.operation == RT_DB_DISCOVERY_TOKEN_SUBTRACT ||
	    node.operation == RT_DB_DISCOVERY_TOKEN_XOR) {
	    if (stack.size() < 2) {
		arcs.clear();
		return;
	    }
	    node.right = stack.back();
	    stack.pop_back();
	    node.left = stack.back();
	    stack.pop_back();
	} else if (node.operation == RT_DB_DISCOVERY_TOKEN_NOT) {
	    if (stack.empty()) {
		arcs.clear();
		return;
	    }
	    node.left = stack.back();
	    stack.pop_back();
	} else {
	    arcs.clear();
	    return;
	}
	nodes.push_back(node);
	stack.push_back(nodes.size() - 1);
    }
    if (stack.size() != 1 || leaf_cursor != leaf_end ||
	arcs.size() != leaf_count) {
	arcs.clear();
	return;
    }

    std::vector<std::pair<size_t, int>> visit;
    visit.emplace_back(stack.back(), OP_UNION);
    while (!visit.empty()) {
	const std::pair<size_t, int> current = visit.back();
	visit.pop_back();
	if (current.first >= nodes.size()) {
	    arcs.clear();
	    return;
	}
	const expression_node &node = nodes[current.first];
	if (node.operation == RT_DB_DISCOVERY_TOKEN_LEAF) {
	    if (node.arc >= arcs.size()) {
		arcs.clear();
		return;
	    }
	    arcs[node.arc].operation =
		rt_db_discovery_public_op(current.second);
	    continue;
	}
	int left_operation = OP_UNION;
	int right_operation = OP_UNION;
	if (node.operation == RT_DB_DISCOVERY_TOKEN_INTERSECT)
	    right_operation = OP_INTERSECT;
	else if (node.operation == RT_DB_DISCOVERY_TOKEN_SUBTRACT)
	    right_operation = OP_SUBTRACT;
	else if (node.operation == RT_DB_DISCOVERY_TOKEN_XOR)
	    right_operation = current.second;
	if (node.right != SIZE_MAX)
	    visit.emplace_back(node.right, right_operation);
	if (node.left != SIZE_MAX)
	    visit.emplace_back(node.left, left_operation);
    }
}


static void
rt_db_discovery_add_primitive_reference(struct db_i *dbip,
	struct directory *parent, const char *name,
	std::vector<rt_db_hierarchy_arc_native> &arcs)
{
    if (!name || !name[0])
	return;
    rt_db_hierarchy_arc_native arc;
    arc.parent = parent;
    arc.child_name = name;
    arc.child = db_lookup(dbip, name, LOOKUP_QUIET);
    arcs.push_back(std::move(arc));
}


static void
rt_db_discovery_decode_primitive(struct db_i *dbip, struct directory *dp,
	std::vector<rt_db_hierarchy_arc_native> &arcs)
{
    if (db_version(dbip) != 5) {
	struct rt_db_internal intern;
	RT_DB_INTERNAL_INIT(&intern);
	if (rt_db_get_internal(&intern, dp, dbip, nullptr) < 0)
	    return;
	if (intern.idb_type == ID_EXTRUDE) {
	    const struct rt_extrude_internal *extrude =
		static_cast<const struct rt_extrude_internal *>(intern.idb_ptr);
	    if (extrude)
		rt_db_discovery_add_primitive_reference(dbip, dp,
		    extrude->sketch_name, arcs);
	} else if (intern.idb_type == ID_REVOLVE) {
	    const struct rt_revolve_internal *revolve =
		static_cast<const struct rt_revolve_internal *>(intern.idb_ptr);
	    if (revolve)
		rt_db_discovery_add_primitive_reference(dbip, dp,
		    bu_vls_cstr(&revolve->sketch_name), arcs);
	} else if (intern.idb_type == ID_DSP) {
	    const struct rt_dsp_internal *dsp =
		static_cast<const struct rt_dsp_internal *>(intern.idb_ptr);
	    if (dsp && dsp->dsp_datasrc == RT_DSP_SRC_OBJ)
		rt_db_discovery_add_primitive_reference(dbip, dp,
		    bu_vls_cstr(&dsp->dsp_name), arcs);
	}
	rt_db_free_internal(&intern);
	return;
    }

    struct bu_external external;
    if (!rt_db_discovery_external_view(dbip, dp, external))
	return;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_external5_to_internal5(&intern, &external, dp->d_namep,
	dbip, nullptr) < 0)
	return;

    if (intern.idb_type == ID_EXTRUDE) {
	const struct rt_extrude_internal *extrude =
	    static_cast<const struct rt_extrude_internal *>(intern.idb_ptr);
	if (extrude)
	    rt_db_discovery_add_primitive_reference(dbip, dp,
		extrude->sketch_name, arcs);
    } else if (intern.idb_type == ID_REVOLVE) {
	const struct rt_revolve_internal *revolve =
	    static_cast<const struct rt_revolve_internal *>(intern.idb_ptr);
	if (revolve)
	    rt_db_discovery_add_primitive_reference(dbip, dp,
		bu_vls_cstr(&revolve->sketch_name), arcs);
    } else if (intern.idb_type == ID_DSP) {
	const struct rt_dsp_internal *dsp =
	    static_cast<const struct rt_dsp_internal *>(intern.idb_ptr);
	if (dsp && dsp->dsp_datasrc == RT_DSP_SRC_OBJ)
	    rt_db_discovery_add_primitive_reference(dbip, dp,
		bu_vls_cstr(&dsp->dsp_name), arcs);
    }
    rt_db_free_internal(&intern);
}


struct rt_db_discovery_worker_state {
    struct db_i *dbip = nullptr;
    const std::vector<rt_db_discovery_task> *tasks = nullptr;
    std::vector<rt_db_discovery_task_result> *results = nullptr;
    std::atomic<size_t> next{0};
};


static void
rt_db_discovery_worker(int UNUSED(cpu), void *data)
{
    rt_db_discovery_worker_state *state =
	static_cast<rt_db_discovery_worker_state *>(data);
    if (!state || !state->dbip || !state->tasks || !state->results)
	return;

    for (;;) {
	const size_t index = state->next.fetch_add(1);
	if (index >= state->tasks->size())
	    return;
	const rt_db_discovery_task &task = (*state->tasks)[index];
	std::vector<rt_db_hierarchy_arc_native> &arcs =
	    (*state->results)[index].arcs;
	if (task.dp && (task.dp->d_flags & RT_DIR_COMB))
	    rt_db_discovery_decode_comb(state->dbip, task.dp, arcs);
	else if (task.dp)
	    rt_db_discovery_decode_primitive(state->dbip, task.dp, arcs);
    }
}


static void
rt_db_discovery_callback(struct db_i *dbip, struct directory *parent,
	struct directory *child, const char *name, db_op_t operation,
	matp_t matrix)
{
    if (!dbip || !dbip->i)
	return;
    for (size_t i = 0;
	i < BU_PTBL_LEN(&dbip->i->dbi_update_nref_clbks); ++i) {
	struct dbi_update_nref_clbk *callback =
	    reinterpret_cast<struct dbi_update_nref_clbk *>(
		BU_PTBL_GET(&dbip->i->dbi_update_nref_clbks, i));
	if (callback && callback->f)
	    callback->f(dbip, parent, child, name, operation, matrix,
		callback->u_data);
    }
}


void
rt_db_discovery_options_init(struct rt_db_discovery_options *options)
{
    if (!options)
	return;
    options->max_workers = 0;
    options->max_inflight_bytes = 0;
}


int
rt_db_discovery_build(struct db_i *dbip,
	const struct rt_db_discovery_options *options,
	struct rt_db_hierarchy **hierarchy,
	struct rt_db_discovery_stats *stats)
{
    if (hierarchy)
	*hierarchy = nullptr;
    if (stats)
	std::memset(stats, 0, sizeof(*stats));
    if (!dbip || !dbip->i || !hierarchy)
	return -1;
    RT_CK_DBI(dbip);

    struct rt_db_discovery_options configured;
    rt_db_discovery_options_init(&configured);
    if (options)
	configured = *options;

    rt_db_discovery_mapped_view mapped(dbip);
    const auto directory_begin = std::chrono::steady_clock::now();
    if (db_dirbuild_without_nref(dbip) < 0)
	return -1;
    if (stats) {
	stats->directory_microseconds =
	    rt_db_discovery_elapsed_us(directory_begin);
	stats->used_mapped_view = mapped.active() ? 1 : 0;
    }

    std::vector<rt_db_discovery_task> tasks;
    struct directory *dp = RT_DIR_NULL;
    uint64_t object_count = 0;
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	object_count++;
	if ((dp->d_flags & RT_DIR_COMB) ||
	    (dp->d_major_type == DB5_MAJORTYPE_BRLCAD &&
	     (dp->d_minor_type == DB5_MINORTYPE_BRLCAD_EXTRUDE ||
	      dp->d_minor_type == DB5_MINORTYPE_BRLCAD_REVOLVE ||
	      dp->d_minor_type == DB5_MINORTYPE_BRLCAD_DSP))) {
	    rt_db_discovery_task task;
	    task.dp = dp;
	    task.serialized_bytes = dp->d_len;
	    tasks.push_back(task);
	}
    } FOR_ALL_DIRECTORY_END;

    size_t workers = configured.max_workers ? configured.max_workers :
	std::min<size_t>(bu_avail_cpus(), 8);
    workers = std::max<size_t>(1, std::min<size_t>(workers, 32));
    workers = std::min<size_t>(workers, std::max<size_t>(1, tasks.size()));
    size_t byte_limit = configured.max_inflight_bytes;
    if (!byte_limit)
	byte_limit = 256ULL * 1024ULL * 1024ULL;
    size_t largest_task = 1;
    for (const rt_db_discovery_task &task : tasks)
	largest_task = std::max(largest_task, task.serialized_bytes);
    /* Each worker owns at most one decoded record.  Bound worker count from
	 * the largest referencing object rather than locking a condition variable
	 * for every one of hundreds of thousands of small combinations. */
    if (largest_task < byte_limit)
	workers = std::min(workers,
	    std::max<size_t>(1, byte_limit / largest_task));
    else
	workers = 1;

    const auto decode_begin = std::chrono::steady_clock::now();
    std::vector<rt_db_discovery_task_result> results(tasks.size());
    rt_db_discovery_worker_state worker_state;
    worker_state.dbip = dbip;
    worker_state.tasks = &tasks;
    worker_state.results = &results;
    if (!tasks.empty())
	bu_parallel(rt_db_discovery_worker, workers, &worker_state);

    rt_db_hierarchy *discovered = new rt_db_hierarchy;
    size_t reference_count = 0;
    for (const rt_db_discovery_task_result &result : results)
	reference_count += result.arcs.size();
    discovered->arcs.reserve(reference_count);
    for (rt_db_discovery_task_result &result : results) {
	for (rt_db_hierarchy_arc_native &arc : result.arcs)
	    discovered->arcs.push_back(std::move(arc));
    }
    if (stats)
	stats->decode_microseconds = rt_db_discovery_elapsed_us(decode_begin);

    const auto publish_begin = std::chrono::steady_clock::now();
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	dp->d_nref = 0;
    } FOR_ALL_DIRECTORY_END;
    rt_db_discovery_callback(dbip, nullptr, nullptr, nullptr,
	DB_OP_UNION, nullptr);
    for (rt_db_hierarchy_arc_native &arc : discovered->arcs) {
	if (arc.child)
	    arc.child->d_nref++;
	rt_db_discovery_callback(dbip, arc.parent, arc.child,
	    arc.child_name.c_str(), arc.operation,
	    arc.matrix_valid ? arc.matrix : nullptr);
    }
    rt_db_discovery_callback(dbip, nullptr, nullptr, nullptr,
	DB_OP_SUBTRACT, nullptr);

    if (stats) {
	stats->publish_microseconds =
	    rt_db_discovery_elapsed_us(publish_begin);
	stats->object_count = object_count;
	stats->referencing_object_count = tasks.size();
	stats->reference_count = discovered->arcs.size();
	stats->worker_count = static_cast<uint32_t>(workers);
    }
    *hierarchy = discovered;
    return 0;
}


void
rt_db_hierarchy_destroy(struct rt_db_hierarchy *hierarchy)
{
    delete hierarchy;
}


size_t
rt_db_hierarchy_arc_count(const struct rt_db_hierarchy *hierarchy)
{
    return hierarchy ? hierarchy->arcs.size() : 0;
}


int
rt_db_hierarchy_arc_get(const struct rt_db_hierarchy *hierarchy,
	size_t index, struct rt_db_hierarchy_arc *arc)
{
    if (!hierarchy || !arc || index >= hierarchy->arcs.size())
	return 0;
    const rt_db_hierarchy_arc_native &native = hierarchy->arcs[index];
    arc->parent = native.parent;
    arc->child = native.child;
    arc->child_name = native.child_name.c_str();
    arc->operation = native.operation;
    arc->matrix_valid = native.matrix_valid ? 1 : 0;
    MAT_COPY(arc->matrix, native.matrix);
    return 1;
}
