# Reviewed installed libBObol header tiers.
#
# Stable headers are included by BObol.h and form the supported host-facing
# source/API contract.  Advanced headers are installed for drawing-aware
# BRL-CAD integrations, but their C++ ABI remains experimental for SOVERSION 1.

# This smaller list is the binary compatibility boundary enforced by the
# strict Doxygen and exported-symbol gates.  The headers are parsed as C, so
# C++ convenience wrappers in mixed-language headers remain source-supported
# rather than becoming part of the SOVERSION promise.
set(LIBBOBOL_STABLE_C_ABI_HEADERS
  BDefines.h
  BDisplayEndpoint.h
  BDisplaySession.h
  BHostFactory.h
  BInput.h
)

set(LIBBOBOL_STABLE_HEADERS
  BDatabaseSource.h
  BDefines.h
  BDisplayEndpoint.h
  BDisplaySession.h
  BExportAction.h
  BFramebuffer.h
  BHeadlessWindowHost.h
  BHostFactory.h
  BInit.h
  BInput.h
  BMeasureAction.h
  BPickDetail.h
  BRtRender.h
  BSceneController.h
  BSceneGroup.h
  BSnapAction.h
  BSourceRealization.h
  BSourceMeshRequest.h
  BViewController.h
  BViewQuery.h
  BViewStore.h
  BWindowHost.h
)

set(LIBBOBOL_ADVANCED_HEADERS
  BADC.h
  BAxes.h
  BDrawCache.h
  BEditPreview.h
  BEvaluatedPoints.h
  BGrid.h
  BHUDLabelOverlay.h
  BImagePlane.h
  BImageSource.h
  BLineLayerOverlay.h
  BLodIdentifiers.h
  BLodMeshShape.h
  BLodRealization.h
  BLodService.h
  BLodUpdateAction.h
  BMaterialObject.h
  BMeshLodCache.h
  BMeshLodSubmitAction.h
  BMeshResidencyAction.h
  BMeshShape.h
  BPerformance.h
  BRealizeAction.h
  BViewAttachment.h
  BViewLod.h
  BViewportImage.h
  BVListShape.h
)

set(LIBBOBOL_INSTALLED_HEADERS
  ${LIBBOBOL_STABLE_HEADERS}
  ${LIBBOBOL_ADVANCED_HEADERS}
)
