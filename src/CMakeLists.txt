add_library(MHtestFunctions mhtestfunctions.cpp)

add_executable(NelderMeadMethod demoNelderMead.cxx)
add_executable(LineSearch demoLineSearch.cxx)
target_link_libraries(NelderMeadMethod MHtestFunctions)
target_link_libraries(LineSearch MHtestFunctions)

install(TARGETS NelderMeadMethod DESTINATION bin)
install(TARGETS MHtestFunctions DESTINATION lib)
