cmake_minimum_required(VERSION 3.1.0)

project(gencut)

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
	"app/GeneMatrix.cpp"
	"app/GeneSignal.cpp"
	"app/GeneSignalFactory.cpp"
	"app/GeneSignalList.cpp"
	"app/GeneSignalMatrix.cpp"
	"app/GenomeCutter.cpp"
	"app/GenomeModel.cpp"
)

add_executable(gencut ${SOURCES})
target_link_libraries(gencut rapidcsv)

