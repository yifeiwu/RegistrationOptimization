CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(MUIGLSPreOp)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
        INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build MUIGLSPreOp without VTK.  Please set VTK_BUILD_DIR.")
ENDIF(VTK_FOUND)

FIND_PACKAGE(VTKMY)
IF(VTKMY_FOUND)
        INCLUDE(${VTKMY_USE_FILE})
ELSE(VTKMY_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build MUIGLSPreOp without without vtkMy. Please set VTKMY_BUILD_DIR.")
ENDIF(VTKMY_FOUND)

# link vtkMy-build/bin
LINK_DIRECTORIES(${VTKMY_DIR}/bin)

LINK_DIRECTORIES(/usr/local/lib)
INCLUDE(/usr/local/lib)
Include_directories(/usr/local/include/eigen3)
INCLUDE_DIRECTORIES(/usr/local/include/eigen3/Eigen)
Include_directories(/home/caleb/src/eigen/test)
INCLUDE_DIRECTORIES (/usr/local/include)
INCLUDE_DIRECTORIES(/home/caleb/src/kdtree-0.5.6)



ADD_EXECUTABLE(IntraCalc IntraCalc.cxx)
TARGET_LINK_LIBRARIES(IntraCalc
     vtkCommon
     vtkGraphics
     vtkHybrid
     vtkRendering
     vtkIO  
     /usr/local/lib/liblmmin.a
     /home/caleb/src/kdtree-0.5.6/libkdtree.a
)
