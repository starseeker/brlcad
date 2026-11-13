/* T R A N S A C T I O N _ F A U L T _ T E S T _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_TESTS_TRANSACTION_FAULT_TEST_PRIVATE_H
#define LIBBOBOL_TESTS_TRANSACTION_FAULT_TEST_PRIVATE_H

#include "../transaction_fault_private.h"

#include "bu/env.h"

#include <cstdlib>
#include <string>

class ScopedTransactionFault {
public:
    explicit ScopedTransactionFault(BObolTransactionFaultPoint point)
    {
	const char *configured = std::getenv(BObolTransactionFaultEnvironment);
	if (configured) {
	    hadPrevious = true;
	    previous = configured;
	}
	bu_setenv(BObolTransactionFaultEnvironment,
	    bobol_transaction_fault_name(point), 1);
    }

    ~ScopedTransactionFault()
    {
	bu_setenv(BObolTransactionFaultEnvironment,
	    hadPrevious ? previous.c_str() : "", 1);
    }

    ScopedTransactionFault(const ScopedTransactionFault &) = delete;
    ScopedTransactionFault &operator=(const ScopedTransactionFault &) = delete;

private:
    bool hadPrevious = false;
    std::string previous;
};

#endif /* LIBBOBOL_TESTS_TRANSACTION_FAULT_TEST_PRIVATE_H */
