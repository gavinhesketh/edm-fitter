cmake_minimum_required(VERSION 3.10)

project(EDMFitter)

# compile the blinding code
add_library(Blinders SHARED Blinders.cc)
target_include_directories(Blinders PUBLIC ${CMAKE_INSTALL_PREFIX}/include)

# and the main executable
include_directories(${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES})
link_directories(${CMAKE_INSTALL_PREFIX}/lib)


set(EXE fitEDM_with_Blinding)
add_executable(${EXE} fitEDM_with_Blinding.C)
#set_target_properties(${EXE} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
target_link_libraries(${EXE} Blinders Random)

#ssl crypto

install(TARGETS Blinders DESTINATION lib)
install(FILES Blinders.hh DESTINATION include)
install(TARGETS ${EXE} DESTINATION bin)
