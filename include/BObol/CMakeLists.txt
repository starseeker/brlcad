set(
  BOBOL_HEADERS
  BADC.h
  BAxes.h
  BDatabaseSource.h
  BDisplaySession.h
  BDefines.h
  BDisplayEndpoint.h
  BDrawCache.h
  BEditPreview.h
  BEvaluatedPoints.h
  BEvaluatedWire.h
  BExportAction.h
  BFramebuffer.h
  BGrid.h
  BHeadlessWindowHost.h
  BHostFactory.h
  BHUDLabelOverlay.h
  BImagePlane.h
  BImageSource.h
  BInit.h
  BInput.h
  BLineLayerOverlay.h
  BLodRealization.h
  BLodService.h
  BLodMeshShape.h
  BLodUpdateAction.h
  BMaterialObject.h
  BMeasureAction.h
  BMeshLodCache.h
  BMeshLodSubmitAction.h
  BMeshResidencyAction.h
  BMeshShape.h
  BPickDetail.h
  BPerformance.h
  BRealizeAction.h
  BRtRender.h
  BSceneController.h
  BSceneGroup.h
  BSnapAction.h
  BSourceMeshRequest.h
  BViewAttachment.h
  BViewController.h
  BViewLod.h
  BViewQuery.h
  BViewStore.h
  BViewportImage.h
  BVListShape.h
  BWindowHost.h
)
brlcad_manage_files(BOBOL_HEADERS ${INCLUDE_DIR}/brlcad/BObol REQUIRED libBObol)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
