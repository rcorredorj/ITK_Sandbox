#-----------------------------------------------------------------------------
# decide what to compile

OPTION(COMPILE_N4BiasCorrectionApp "Compile N4BiasFieldCorrectionApp" ON)
IF(COMPILE_N4BiasCorrectionApp)
  SUBDIRS(N4BiasFieldCorrectionApp)
ENDIF(COMPILE_N4BiasCorrectionApp)
