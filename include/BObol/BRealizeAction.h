/*              B R E A L I Z E A C T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BRealizeAction.h */

#ifndef BOBOL_BREALIZEACTION_H
#define BOBOL_BREALIZEACTION_H

#include "BObol/BDefines.h"

#include <Inventor/SbString.h>
#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoSubAction.h>
#include <memory>
#include <vector>

class SoBRLDatabaseSource;
class BObolSceneController;
struct BObolDatabaseSourceRealizationCache;
struct BObolSourceRealizationEffects;

class BOBOL_EXPORT BObolRealizationRepository {
public:
    BObolRealizationRepository(void);
    ~BObolRealizationRepository(void);
    BObolRealizationRepository(
	const BObolRealizationRepository &) = delete;
    BObolRealizationRepository &operator=(
	const BObolRealizationRepository &) = delete;

    void clear(void);
    void invalidateObject(const char *name);
    void renameObject(const char *oldName, const char *newName);
    void invalidateViewVariants(void);
    void seedSource(SoBRLDatabaseSource *source);
    void releaseSource(SoBRLDatabaseSource *source);

private:
    friend class SoBRLRealizeAction;
    friend class BObolSceneController;
    struct SourceMembership {
	const BObolSceneController *controller;
	SoBRLDatabaseSource *source;
	bool present;
    };
    class SourceUpdate {
    public:
	~SourceUpdate();
	void commit() noexcept;
    private:
	friend class BObolRealizationRepository;
	struct Impl;
	explicit SourceUpdate(std::unique_ptr<Impl> prepared);
	std::unique_ptr<Impl> d;
    };
    class ObjectRename {
    public:
	~ObjectRename();
	void commit() noexcept;
    private:
	friend class BObolRealizationRepository;
	struct Impl;
	explicit ObjectRename(std::unique_ptr<Impl> prepared);
	std::unique_ptr<Impl> d;
    };
    std::unique_ptr<SourceUpdate> prepareSourceRelease(const std::vector<SoBRLDatabaseSource *> &sources);
    std::unique_ptr<SourceUpdate> prepareSourceSeed(const std::vector<SoBRLDatabaseSource *> &sources);
    std::unique_ptr<SourceUpdate> prepareSourceMembership(
	const std::vector<SourceMembership> &memberships,
	const std::vector<SoBRLDatabaseSource *> &seeded);
    std::unique_ptr<ObjectRename> prepareObjectRename(
	const char *oldName, const char *newName);
    void attachController() noexcept;
    void detachController(const BObolSceneController *controller) noexcept;
    bool hasSourceOwner(const BObolSceneController *controller, const SoBRLDatabaseSource *source) const noexcept;
    bool acceptsSource(const SoBRLDatabaseSource *source) const noexcept;
    struct Residency;
    std::unique_ptr<BObolDatabaseSourceRealizationCache> cache;
    std::unique_ptr<Residency> residency;
};

class BOBOL_EXPORT SoBRLRealizeAction : public SoAction {
    typedef SoAction inherited;

    SO_ACTION_HEADER(SoBRLRealizeAction);

public:
    SoBRLRealizeAction(void);
    virtual ~SoBRLRealizeAction(void);
    static void initClass(void);

    /* Progress belongs to the current/last apply. Completed source counts and
     * diagnostics publish before source observers and survive traversal errors. */
    unsigned int getVisitedSourceCount(void) const;
    unsigned int getRealizedSourceCount(void) const;
    unsigned int getFailedSourceCount(void) const;
    const SbString &getDiagnostics(void) const;
    void setRetainRealizationCache(SbBool retain);
    void clearRealizationCache(void);
    void invalidateRealizationObject(const char *name);
    /* The repository is borrowed and must outlive this action. */
    void setRealizationRepository(BObolRealizationRepository *repository);

protected:
    virtual void beginTraversal(SoNode *node);

private:
    friend class BObolSceneController;
    struct SourcePublication;
    static void nodeAction(SoAction *action, SoNode *node);
    static void databaseSourceAction(SoAction *action, SoNode *node);
    void stopSceneTraversal() noexcept;

    /* Borrowed only while the scene applies this action. */
    SoBRLRealizeAction *enclosingSceneAction = nullptr;
    BObolSourceRealizationEffects *publicationEffects = nullptr;

    unsigned int visitedSourceCount;
    unsigned int realizedSourceCount;
    unsigned int failedSourceCount;
    SbString diagnostics;
    BObolDatabaseSourceRealizationCache *realizationCache;
    BObolRealizationRepository *realizationRepository;
    SbBool ownsRealizationRepository;
    SbBool seedingCache;
    SbBool retainRealizationCache;
};

#endif /* BOBOL_BREALIZEACTION_H */
