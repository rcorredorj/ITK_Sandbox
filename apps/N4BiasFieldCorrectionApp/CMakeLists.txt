
SET(N4BiasFieldCorrectionApp_SRCS
  main.cpp
)

ADD_EXECUTABLE(N4BiasFieldCorrectionApp ${N4BiasFieldCorrectionApp_SRCS})

TARGET_LINK_LIBRARIES(N4BiasFieldCorrectionApp ${ITK_LIBRARIES})
