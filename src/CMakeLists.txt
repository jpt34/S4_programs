cmake_minimum_required(VERSION 3.5)

project(S4_programs LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

set( CMAKE_CXX_FLAGS "-std=c++1z -mavx -fPIC")

if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message("Release build")
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBS OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

add_executable(S4_programs main.cpp Mlogit.cpp bgen/IndexQuery.cpp bgen/MissingValue.cpp bgen/View.cpp bgen/bgen.cpp bgen/ImputeclassDetailed.cpp bgen/db/Connection.cpp bgen/db/Error.cpp bgen/db/SQLStatement.cpp bgen/db/SQLite3Connection.cpp bgen/db/SQLite3Statement.cpp routines.cpp)

target_link_libraries(S4_programs -lboost_program_options -lboost_iostreams -lz -lzstd -lboost_system -lboost_filesystem -lboost_thread -pthread /home/jpt34/ExternalLibrary/libboost.a -lrt -lsqlite3 -dynamic -ldl)
