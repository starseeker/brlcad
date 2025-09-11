/*                T E S T _ D B _ T H R E A D I N G . C
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
/** @file test_db_threading.c
 *
 * Comprehensive multithreaded stress test for BRL-CAD database I/O operations.
 * Tests that the C++17-based locking mechanism correctly protects database
 * operations under various concurrent access patterns.
 *
 */

#include "common.h"

#include <string.h>
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#include "bio.h"

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/snooze.h"
#include "bu/str.h"
#include "bu/time.h"
#include "vmath.h"
#include "raytrace.h"
#include "wdb.h"

/* Test configuration */
#define MAX_THREADS 8
#define OBJECTS_PER_THREAD 50
#define TEST_ITERATIONS 5
#define MAX_OBJECT_NAME 64

/* Thread test data structure */
struct thread_data {
    struct db_i *dbip;          /* Database instance */
    struct rt_wdb *wdbp;        /* Write database instance */
    int thread_id;              /* Unique thread identifier */
    int test_mode;              /* Test mode (1=read, 2=write, 3=mixed) */
    int object_count;           /* Number of objects to process */
    int success_count;          /* Number of successful operations */
    int error_count;            /* Number of failed operations */
    char test_name[64];         /* Test description */
};

/* Global test results */
static int g_total_success = 0;
static int g_total_errors = 0;
static int g_test_failures = 0;

/**
 * Create a sphere object in the database
 */
static int
create_sphere_object(struct rt_wdb *wdbp, const char *name, point_t center, fastf_t radius)
{
    if (!wdbp || !name) {
        return -1;
    }
    return mk_sph(wdbp, name, center, radius);
}

/**
 * Create an ARB8 (box) object in the database
 */
static int
create_arb8_object(struct rt_wdb *wdbp, const char *name, point_t *vertices)
{
    if (!wdbp || !name || !vertices) {
        return -1;
    }
    return mk_arb8(wdbp, name, (const fastf_t *)vertices);
}

/**
 * Read and validate an object from the database
 */
static int
read_and_validate_object(struct db_i *dbip, const char *name)
{
    struct directory *dp;
    struct rt_db_internal intern;
    
    if (!dbip || !name) {
        return -1;
    }
    
    dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
        return -1; /* Object not found */
    }
    
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, bn_mat_identity, &rt_uniresource) < 0) {
        return -1;
    }
    
    /* Validate that we got a valid object */
    if (intern.idb_ptr == NULL) {
        rt_db_free_internal(&intern);
        return -1;
    }
    
    rt_db_free_internal(&intern);
    return 0;
}

/**
 * Delete an object from the database
 */
static int
delete_object(struct db_i *dbip, const char *name)
{
    struct directory *dp;
    
    if (!dbip || !name) {
        return -1;
    }
    
    dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
        return -1; /* Object not found */
    }
    
    return db_delete(dbip, dp);
}

/**
 * Thread worker function - performs concurrent read operations
 */
static void
reader_thread_worker(int cpu_num, void *arg)
{
    struct thread_data *td = (struct thread_data *)arg;
    char obj_name[MAX_OBJECT_NAME];
    int i;
    
    bu_log("Reader thread %d starting (%s)\n", td->thread_id, td->test_name);
    
    for (i = 0; i < td->object_count; i++) {
        /* Try to read various objects created by other threads */
        snprintf(obj_name, sizeof(obj_name), "sphere_%d_%d", (i % MAX_THREADS), i);
        
        if (read_and_validate_object(td->dbip, obj_name) == 0) {
            td->success_count++;
        } else {
            td->error_count++;
        }
        
        /* Also try reading ARB8 objects */
        snprintf(obj_name, sizeof(obj_name), "arb8_%d_%d", (i % MAX_THREADS), i);
        if (read_and_validate_object(td->dbip, obj_name) == 0) {
            td->success_count++;
        } else {
            td->error_count++;
        }
        
        /* Small delay to allow other threads to work */
        bu_snooze(1000); /* 1ms */
    }
    
    bu_log("Reader thread %d completed: %d successes, %d errors\n", 
           td->thread_id, td->success_count, td->error_count);
}

/**
 * Thread worker function - performs concurrent write operations
 */
static void
writer_thread_worker(int cpu_num, void *arg)
{
    struct thread_data *td = (struct thread_data *)arg;
    char obj_name[MAX_OBJECT_NAME];
    point_t center;
    point_t vertices[8];
    fastf_t radius;
    int i, j;
    
    bu_log("Writer thread %d starting (%s)\n", td->thread_id, td->test_name);
    
    for (i = 0; i < td->object_count; i++) {
        /* Create sphere objects */
        snprintf(obj_name, sizeof(obj_name), "sphere_%d_%d", td->thread_id, i);
        VSET(center, td->thread_id * 10.0, i * 5.0, 0.0);
        radius = 1.0 + (td->thread_id * 0.1) + (i * 0.01);
        
        if (create_sphere_object(td->wdbp, obj_name, center, radius) == 0) {
            td->success_count++;
        } else {
            td->error_count++;
        }
        
        /* Create ARB8 objects */
        snprintf(obj_name, sizeof(obj_name), "arb8_%d_%d", td->thread_id, i);
        
        /* Generate vertices for a box */
        for (j = 0; j < 8; j++) {
            VSET(vertices[j], 
                 td->thread_id * 15.0 + ((j & 1) ? 2.0 : 0.0),
                 i * 8.0 + ((j & 2) ? 2.0 : 0.0),
                 ((j & 4) ? 2.0 : 0.0));
        }
        
        if (create_arb8_object(td->wdbp, obj_name, vertices) == 0) {
            td->success_count++;
        } else {
            td->error_count++;
        }
        
        /* Small delay to allow other threads to work */
        bu_snooze(2000); /* 2ms */
    }
    
    bu_log("Writer thread %d completed: %d successes, %d errors\n", 
           td->thread_id, td->success_count, td->error_count);
}

/**
 * Thread worker function - performs mixed read/write/delete operations
 */
static void
mixed_thread_worker(int cpu_num, void *arg)
{
    struct thread_data *td = (struct thread_data *)arg;
    char obj_name[MAX_OBJECT_NAME];
    char temp_name[MAX_OBJECT_NAME];
    point_t center;
    fastf_t radius;
    int i;
    
    bu_log("Mixed thread %d starting (%s)\n", td->thread_id, td->test_name);
    
    for (i = 0; i < td->object_count; i++) {
        /* Create a temporary object */
        snprintf(temp_name, sizeof(temp_name), "temp_%d_%d", td->thread_id, i);
        VSET(center, td->thread_id * 20.0, i * 10.0, 5.0);
        radius = 0.5 + (td->thread_id * 0.05) + (i * 0.005);
        
        if (create_sphere_object(td->wdbp, temp_name, center, radius) == 0) {
            td->success_count++;
            
            /* Try to read it back immediately */
            if (read_and_validate_object(td->dbip, temp_name) == 0) {
                td->success_count++;
                
                /* Try to delete it */
                if (delete_object(td->dbip, temp_name) == 0) {
                    td->success_count++;
                } else {
                    td->error_count++;
                }
            } else {
                td->error_count++;
            }
        } else {
            td->error_count++;
        }
        
        /* Try to read objects created by other threads */
        snprintf(obj_name, sizeof(obj_name), "sphere_%d_%d", (i % MAX_THREADS), i / 2);
        if (read_and_validate_object(td->dbip, obj_name) == 0) {
            td->success_count++;
        } else {
            td->error_count++;
        }
        
        /* Small delay */
        bu_snooze(1500); /* 1.5ms */
    }
    
    bu_log("Mixed thread %d completed: %d successes, %d errors\n", 
           td->thread_id, td->success_count, td->error_count);
}

/**
 * Run a multithreaded test with specified configuration
 */
static int
run_threaded_test(struct db_i *dbip, struct rt_wdb *wdbp, 
                  const char *test_name, int test_mode, 
                  int num_threads, int objects_per_thread)
{
    struct thread_data thread_data[MAX_THREADS];
    void (*worker_func)(int, void *);
    int i;
    int total_success = 0;
    int total_errors = 0;
    
    if (num_threads > MAX_THREADS) {
        bu_log("ERROR: Too many threads requested (%d > %d)\n", num_threads, MAX_THREADS);
        return -1;
    }
    
    bu_log("\n=== Starting Test: %s ===\n", test_name);
    bu_log("Threads: %d, Objects per thread: %d\n", num_threads, objects_per_thread);
    
    /* Choose worker function based on test mode */
    switch (test_mode) {
        case 1: /* Read-only test */
            worker_func = reader_thread_worker;
            break;
        case 2: /* Write-only test */
            worker_func = writer_thread_worker;
            break;
        case 3: /* Mixed read/write/delete test */
            worker_func = mixed_thread_worker;
            break;
        default:
            bu_log("ERROR: Invalid test mode: %d\n", test_mode);
            return -1;
    }
    
    /* Initialize thread data */
    for (i = 0; i < num_threads; i++) {
        thread_data[i].dbip = dbip;
        thread_data[i].wdbp = wdbp;
        thread_data[i].thread_id = i;
        thread_data[i].test_mode = test_mode;
        thread_data[i].object_count = objects_per_thread;
        thread_data[i].success_count = 0;
        thread_data[i].error_count = 0;
        bu_strlcpy(thread_data[i].test_name, test_name, sizeof(thread_data[i].test_name));
    }
    
    /* Run the threads */
    bu_parallel(worker_func, num_threads, (void *)thread_data);
    
    /* Collect results */
    for (i = 0; i < num_threads; i++) {
        total_success += thread_data[i].success_count;
        total_errors += thread_data[i].error_count;
    }
    
    bu_log("=== Test Results: %s ===\n", test_name);
    bu_log("Total successes: %d\n", total_success);
    bu_log("Total errors: %d\n", total_errors);
    
    /* Update global counters */
    g_total_success += total_success;
    g_total_errors += total_errors;
    
    if (total_errors > 0) {
        g_test_failures++;
        bu_log("FAILED: %s had %d errors\n", test_name, total_errors);
        return -1;
    } else {
        bu_log("PASSED: %s completed successfully\n", test_name);
        return 0;
    }
}

/**
 * Create initial test objects for read tests
 */
static int
create_initial_test_objects(struct rt_wdb *wdbp)
{
    char obj_name[MAX_OBJECT_NAME];
    point_t center;
    point_t vertices[8];
    fastf_t radius;
    int i, j, k;
    int success_count = 0;
    
    bu_log("Creating initial test objects...\n");
    
    /* Create sphere objects for each thread to read */
    for (i = 0; i < MAX_THREADS; i++) {
        for (j = 0; j < OBJECTS_PER_THREAD; j++) {
            snprintf(obj_name, sizeof(obj_name), "sphere_%d_%d", i, j);
            VSET(center, i * 100.0, j * 50.0, 0.0);
            radius = 5.0 + (i * 0.5) + (j * 0.1);
            
            if (create_sphere_object(wdbp, obj_name, center, radius) == 0) {
                success_count++;
            }
            
            /* Create ARB8 objects too */
            snprintf(obj_name, sizeof(obj_name), "arb8_%d_%d", i, j);
            for (k = 0; k < 8; k++) {
                VSET(vertices[k], 
                     i * 105.0 + ((k & 1) ? 10.0 : 0.0),
                     j * 55.0 + ((k & 2) ? 10.0 : 0.0),
                     ((k & 4) ? 10.0 : 0.0));
            }
            
            if (create_arb8_object(wdbp, obj_name, vertices) == 0) {
                success_count++;
            }
        }
    }
    
    bu_log("Created %d initial test objects\n", success_count);
    return success_count;
}

/**
 * Main test function
 */
int
main(int argc, char **argv)
{
    struct db_i *dbip;
    struct rt_wdb *wdbp;
    char *test_db_file = "test_threading.g";
    int result = 0;
    
    bu_setprogname(argv[0]);
    
    /* Clean up any existing test database */
    if (bu_file_exists(test_db_file, NULL)) {
        bu_file_delete(test_db_file);
    }
    
    bu_log("=== BRL-CAD Database Threading Stress Test ===\n");
    bu_log("Testing C++17-based read/write locking implementation\n");
    bu_log("Database file: %s\n", test_db_file);
    
    /* Create a new database */
    wdbp = wdb_fopen(test_db_file);
    if (!wdbp) {
        bu_log("FATAL: Cannot create test database: %s\n", test_db_file);
        return 1;
    }
    
    dbip = wdbp->dbip;
    
    /* Set up database */
    mk_id(wdbp, "Threading Stress Test Database");
    
    /* Create initial objects for read tests */
    create_initial_test_objects(wdbp);
    
    /* Test 1: Multiple concurrent readers */
    if (run_threaded_test(dbip, wdbp, "Concurrent Readers", 1, 4, OBJECTS_PER_THREAD) != 0) {
        result = 1;
    }
    
    /* Test 2: Multiple concurrent writers */
    if (run_threaded_test(dbip, wdbp, "Concurrent Writers", 2, 3, OBJECTS_PER_THREAD / 2) != 0) {
        result = 1;
    }
    
    /* Test 3: Mixed concurrent operations */
    if (run_threaded_test(dbip, wdbp, "Mixed Read/Write/Delete", 3, 4, OBJECTS_PER_THREAD / 3) != 0) {
        result = 1;
    }
    
    /* Test 4: Heavy load test */
    if (run_threaded_test(dbip, wdbp, "Heavy Load Test", 3, 6, OBJECTS_PER_THREAD / 4) != 0) {
        result = 1;
    }
    
    /* Clean up */
    wdb_close(wdbp);
    
    /* Print final results */
    bu_log("\n=== FINAL TEST RESULTS ===\n");
    bu_log("Total successful operations: %d\n", g_total_success);
    bu_log("Total failed operations: %d\n", g_total_errors);
    bu_log("Failed test suites: %d\n", g_test_failures);
    
    if (result == 0 && g_total_errors == 0 && g_test_failures == 0) {
        bu_log("✅ ALL TESTS PASSED: Database threading is working correctly!\n");
    } else {
        bu_log("❌ TESTS FAILED: %d test suites failed with %d total errors\n", 
               g_test_failures, g_total_errors);
    }
    
    /* Clean up test database */
    bu_file_delete(test_db_file);
    
    return result;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */