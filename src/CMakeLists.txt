cmake_minimum_required(VERSION 3.1.0)

project(gencut)

add_compile_definitions(SKIP_R_API)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(INCLUDE_PATHS ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${INCLUDE_PATHS})

# need this to include the auto generated ui_mainwindow.h
# file from the qt5_wrap_ui call below.
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

include(FetchContent)
FetchContent_Declare(
  rapidcsv
  GIT_REPOSITORY "https://github.com/d99kris/rapidcsv.git"
  GIT_TAG        "v8.50"
)
FetchContent_MakeAvailable(rapidcsv)

set(SOURCES 
	"gencut.cpp"
	"app/Application.cpp"
	"app/Application.h"
	"app/GeneMatrix.cpp"
	"app/GeneMatrix.h"
	"app/GeneSignal.cpp"
	"app/GeneSignal.h"
	"app/GeneSignalFactory.cpp"
	"app/GeneSignalFactory.h"
	"app/GeneSignalList.cpp"
	"app/GeneSignalList.h"
	"app/GeneSignalMatrix.cpp"
	"app/GeneSignalMatrix.h"
	"app/GenomeCutter.cpp"
	"app/GenomeCutter.h"
	"app/GenomeModel.cpp"
	"app/GenomeModel.h"
	"app/StringConstants.cpp"
	"app/StringConstants.h"
)

add_executable(gencut ${SOURCES})
target_link_libraries(gencut rapidcsv)

