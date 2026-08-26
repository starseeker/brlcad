/*       B P R E S E N T A T I O N _ P R E P A R A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BPresentationPreparation.h */

#ifndef BOBOL_BPRESENTATIONPREPARATION_H
#define BOBOL_BPRESENTATIONPREPARATION_H

/** Aggregate result of one exact-target retained preparation observation. */
enum BObolCadPreparationProgress {
    BOBOL_CAD_PREPARATION_NONE = 0,
    BOBOL_CAD_PREPARATION_STARTED,
    BOBOL_CAD_PREPARATION_ADVANCED,
    BOBOL_CAD_PREPARATION_COMPLETED,
    BOBOL_CAD_PREPARATION_CONSTRAINED,
    BOBOL_CAD_PREPARATION_FAILED
};

#endif /* BOBOL_BPRESENTATIONPREPARATION_H */
