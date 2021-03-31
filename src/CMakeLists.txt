add_library(MHtestFunctions mhtestfunctions.cpp)

add_executable(NelderMeadMethod demoNelderMead.cxx)
target_link_libraries(NelderMeadMethod MHtestFunctions)

install(TARGETS NelderMeadMethod DESTINATION bin)
install(TARGETS MHtestFunctions DESTINATION lib)
