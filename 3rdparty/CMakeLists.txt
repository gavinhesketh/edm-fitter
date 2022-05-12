cmake_minimum_required(VERSION 3.1)

project(EDMFitter)

include(ExternalProject)

#set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_CURRENT_BINARY_DIR})

set(RL_PREFIX RandomLib-1.10)

# RandomLib, needed for blinding:
ExternalProject_add(${RL_PREFIX}
  PREFIX ${RL_PREFIX}
  URL    https://downloads.sourceforge.net/project/randomlib/distrib/RandomLib-1.10.tar.gz
  UPDATE_COMMAND    ""
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
)
