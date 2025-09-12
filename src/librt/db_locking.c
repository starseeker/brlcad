/*                   D B _ L O C K I N G . C
 * BRL-CAD
 *
 * Copyright (c) 2025 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file librt/db_locking.c
 *
 * Pure C intra-process read/write locking for database instances using libbu semaphores.
 * Provides thread-safe read/write access to database I/O operations.
 *
 */

#include "common.h"

#include "bu/parallel.h"
#include "rt/db_instance.h"
#include "librt_private.h"

__BEGIN_DECLS

/**
 * Internal locking structure using libbu semaphores for read/write locking.
 * Implements a classic readers-writer lock using two semaphores:
 * - reader_count_sem: protects the reader count variable
 * - write_sem: provides exclusive access for writers
 */
struct db_lock_internal {
    uint32_t magic;
    int reader_count;           /* Number of active readers */
    int reader_count_sem;       /* Semaphore to protect reader_count */
    int write_sem;              /* Semaphore for exclusive writer access */
};
#define DB_LOCK_MAGIC 0x646C6B00  /* "dlk\0" */

/**
 * Initialize database locking for a database instance.
 * Creates the internal locking structure using libbu semaphores.
 *
 * @param dbip Database instance pointer
 * @return 0 on success, -1 on failure
 */
int
db_lock_init(struct db_i *dbip)
{
    struct db_lock_internal *lock_data;

    if (!dbip || !dbip->i) {
        return -1;
    }

    /* Allocate the locking structure */
    lock_data = (struct db_lock_internal *)bu_calloc(1, sizeof(struct db_lock_internal), "db_lock_internal");
    if (!lock_data) {
        return -1;
    }

    lock_data->magic = DB_LOCK_MAGIC;
    lock_data->reader_count = 0;
    
    /* Register semaphores for this database instance */
    lock_data->reader_count_sem = bu_semaphore_register("db_reader_count");
    lock_data->write_sem = bu_semaphore_register("db_write");
    
    if (lock_data->reader_count_sem < 0 || lock_data->write_sem < 0) {
        bu_free(lock_data, "db_lock_internal");
        return -1;
    }

    /* Initialize the semaphores */
    bu_semaphore_init(lock_data->reader_count_sem + 1);
    bu_semaphore_init(lock_data->write_sem + 1);

    dbip->i->lock_data = (void *)lock_data;
    return 0;
}

/**
 * Destroy database locking for a database instance.
 * Cleans up the internal locking structure.
 *
 * @param dbip Database instance pointer
 */
void
db_lock_destroy(struct db_i *dbip)
{
    struct db_lock_internal *lock_data;

    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    lock_data = (struct db_lock_internal *)dbip->i->lock_data;
    BU_CKMAG(lock_data, DB_LOCK_MAGIC, "db_lock_internal");

    /* Note: bu_semaphore_free() will be called by libbu cleanup, 
     * so we don't need to explicitly free individual semaphores */
    
    bu_free(lock_data, "db_lock_internal");
    dbip->i->lock_data = NULL;
}

/**
 * Acquire a read lock on the database.
 * Multiple threads can hold read locks simultaneously.
 * Implements readers-writer lock algorithm:
 * 1. Acquire reader count mutex
 * 2. Increment reader count
 * 3. If first reader, acquire write mutex (blocks writers)
 * 4. Release reader count mutex
 *
 * @param dbip Database instance pointer
 */
void
db_acquire_read_lock(struct db_i *dbip)
{
    struct db_lock_internal *lock_data;

    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    lock_data = (struct db_lock_internal *)dbip->i->lock_data;
    BU_CKMAG(lock_data, DB_LOCK_MAGIC, "db_lock_internal");

    /* Acquire the reader count semaphore */
    bu_semaphore_acquire(lock_data->reader_count_sem);
    
    /* Increment reader count */
    lock_data->reader_count++;
    
    /* If this is the first reader, acquire the write semaphore to block writers */
    if (lock_data->reader_count == 1) {
        bu_semaphore_acquire(lock_data->write_sem);
    }
    
    /* Release the reader count semaphore */
    bu_semaphore_release(lock_data->reader_count_sem);
}

/**
 * Acquire a write lock on the database.
 * Only one thread can hold a write lock, and no read locks can be held.
 * Simply acquires the write semaphore for exclusive access.
 *
 * @param dbip Database instance pointer
 */
void
db_acquire_write_lock(struct db_i *dbip)
{
    struct db_lock_internal *lock_data;

    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    lock_data = (struct db_lock_internal *)dbip->i->lock_data;
    BU_CKMAG(lock_data, DB_LOCK_MAGIC, "db_lock_internal");

    /* Acquire exclusive write access */
    bu_semaphore_acquire(lock_data->write_sem);
}

/**
 * Release a lock on the database.
 * Works for both read and write locks by checking reader count.
 * For read locks:
 * 1. Acquire reader count mutex
 * 2. Decrement reader count
 * 3. If last reader, release write mutex (allows writers)
 * 4. Release reader count mutex
 * For write locks:
 * Simply releases the write mutex.
 *
 * @param dbip Database instance pointer
 */
void
db_release_lock(struct db_i *dbip)
{
    struct db_lock_internal *lock_data;

    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    lock_data = (struct db_lock_internal *)dbip->i->lock_data;
    BU_CKMAG(lock_data, DB_LOCK_MAGIC, "db_lock_internal");

    /* Try to acquire reader count mutex to see if this is a read lock */
    bu_semaphore_acquire(lock_data->reader_count_sem);
    
    if (lock_data->reader_count > 0) {
        /* This is a read lock being released */
        lock_data->reader_count--;
        
        /* If this was the last reader, release the write semaphore */
        if (lock_data->reader_count == 0) {
            bu_semaphore_release(lock_data->write_sem);
        }
        
        bu_semaphore_release(lock_data->reader_count_sem);
    } else {
        /* This is a write lock being released */
        bu_semaphore_release(lock_data->reader_count_sem);
        bu_semaphore_release(lock_data->write_sem);
    }
}

__END_DECLS

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */