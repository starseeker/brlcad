/*                   D B _ L O C K I N G . C P P
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
/** @file librt/db_locking.cpp
 *
 * C++17-based intra-process read/write locking for database instances.
 * Provides thread-safe read/write access to database I/O operations.
 *
 */

#include "common.h"

#include <memory>
#include <shared_mutex>

#include "rt/db_instance.h"
#include "librt_private.h"

__BEGIN_DECLS

/**
 * Internal C++ locking structure hidden from C API.
 * Uses C++17 shared_mutex for efficient read/write locking.
 */
struct db_lock_internal {
    std::shared_mutex mutex;
};

/**
 * Initialize database locking for a database instance.
 * Creates the internal C++ locking structure.
 *
 * @param dbip Database instance pointer
 * @return 0 on success, -1 on failure
 */
int
db_lock_init(struct db_i *dbip)
{
    if (!dbip || !dbip->i) {
        return -1;
    }

    try {
        dbip->i->lock_data = new db_lock_internal();
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Destroy database locking for a database instance.
 * Cleans up the internal C++ locking structure.
 *
 * @param dbip Database instance pointer
 */
void
db_lock_destroy(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    try {
        db_lock_internal *lock_data = static_cast<db_lock_internal*>(dbip->i->lock_data);
        delete lock_data;
        dbip->i->lock_data = nullptr;
    } catch (...) {
        // Ignore exceptions during cleanup
    }
}

/**
 * Acquire a read lock on the database.
 * Multiple threads can hold read locks simultaneously.
 *
 * @param dbip Database instance pointer
 */
void
db_acquire_read_lock(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    try {
        db_lock_internal *lock_data = static_cast<db_lock_internal*>(dbip->i->lock_data);
        lock_data->mutex.lock_shared();
    } catch (...) {
        // If locking fails, continue without locking (non-blocking)
    }
}

/**
 * Acquire a write lock on the database.
 * Only one thread can hold a write lock, and no read locks can be held.
 *
 * @param dbip Database instance pointer
 */
void
db_acquire_write_lock(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    try {
        db_lock_internal *lock_data = static_cast<db_lock_internal*>(dbip->i->lock_data);
        lock_data->mutex.lock();
    } catch (...) {
        // If locking fails, continue without locking (non-blocking)
    }
}

/**
 * Release a lock on the database.
 * Works for both read and write locks.
 *
 * @param dbip Database instance pointer
 */
void
db_release_lock(struct db_i *dbip)
{
    if (!dbip || !dbip->i || !dbip->i->lock_data) {
        return;
    }

    try {
        db_lock_internal *lock_data = static_cast<db_lock_internal*>(dbip->i->lock_data);
        lock_data->mutex.unlock();
    } catch (...) {
        // If unlocking fails, continue (non-blocking)
    }
}

__END_DECLS

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */