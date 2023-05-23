#################################################################
#                       ANEUSURFISOLATION                       #
#################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

#---------------------------------
# Project name
#---------------------------------

PROJECT(APP)

# Don't build anything unless the version of CMake is high enough.
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 3.0)

    IF(NOT CMAKE_BUILD_TYPE)
        SET(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release
        RelWithDebInfo MinSizeRel." FORCE)
    ENDIF()
    
    ADD_DEFINITIONS(-Wall )     #all warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")    
    
    IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local"
            CACHE PATH "Default install path" FORCE)
    ENDIF()
        
    #Compiler Options 
    SET(CMAKE_CXX_FLAGS_DEBUG           "-pg -g"        CACHE STRING "Flag to generate debug builds"    FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL      "-Os -DNDEBUG " CACHE STRING "Flag to generate minsize builds"  FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE         "-O3 -DNDEBUG"  CACHE STRING "Flag to generate release builds"  FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -g"        CACHE STRING "Flag to generate release builds with debug information" FORCE)

    IF(UNIX)
        # This includes all libraries (.lib) that Zlib generates.
        SET(ZLIB_LIBRARIES z)

        # MLIB this includes all libraries (.lib)
        SET(MLIB_LIBRARIES m)

        #------------------------------------------------------------
        # Link the executable to the necessary libraries.
        SET(APP_LINK_LIBRARIES
            ${ZLIB_LIBRARIES}
            ${MLIB_LIBRARIES}
        )
        
        # include and library directories
        SET(USERLIB_INCLUDE /usr/include)
        SET(USERLIB_LIBRARIES_DIR /usr/lib)
        INCLUDE_DIRECTORIES(${USERLIB_INCLUDE})
        LINK_DIRECTORIES(${USERLIB_LIBRARIES_DIR})
    ENDIF(UNIX)
    
    # VTK includes and libs directories
    FIND_PACKAGE(VTK)
    IF (NOT VTK_FOUND)
        MESSAGE(FATAL_ERROR "Cannot build without VTK. Please set VTK_DIR.")
        RETURN()
    ENDIF()
    
   
    IF(VTK_FOUND AND VTK_VERSION VERSION_GREATER_EQUAL "8.2.0")
    
        SET(APP_LINK_LIBRARIES
                ${VTK_LIBRARIES}
                ${APP_LINK_LIBRARIES}
            )           
        
        include_directories(utilities)        
        
        # "AneuSurfIsolation.cpp"
        ADD_EXECUTABLE(AneuSurfIsolation AneuSurfIsolation.cpp
                                        ./utilities/aneuFunctions.cpp
                                        ./utilities/aneuConfigFile.cpp
                                        ./utilities/auxAneuSurfIsolation.cpp                                        
                                        )
        TARGET_LINK_LIBRARIES(AneuSurfIsolation ${APP_LINK_LIBRARIES})
            
        INSTALL_TARGETS(/bin AneuSurfIsolation)
        
    ELSEIF(VTK_FOUND AND VTK_VERSION VERSION_LESS "8.2.0")
        MESSAGE("The detected version of VTK is too old. Please, donwload version 8.2.0 or higher")
        RETURN()
    ENDIF(VTK_FOUND AND VTK_VERSION VERSION_GREATER_EQUAL "8.2.0")    

# If the version of CMake was too old, complain and build nothing.
# These should be the last lines in this file.
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 3.0)
  MESSAGE(SEND_ERROR
          "This version of CMake is too old to build this app.  "
          "Please upgrade to CMake 3.0 or greater.")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 3.0)


