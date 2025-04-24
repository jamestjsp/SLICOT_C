# CMakeLists.txt for SLICOT main sources (src)
include(GNUInstallDirs) # Ensure install paths are available
include(${CMAKE_SOURCE_DIR}/cmake/slicot_sources.cmake) # Include the source list file

# Create the slicot library with Fortran sources only
# Use the variable defined in slicot_sources.cmake
add_library(slicot ${SLICOT_CORE_SOURCES})

# Set library properties
set_target_properties(slicot PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    POSITION_INDEPENDENT_CODE ON  # Important for shared libraries
    OUTPUT_NAME "slicot") # Keep the base name simple

# Add include directories if needed (e.g., if there were .h files)
# target_include_directories(slicot PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> # Needed if headers are in src
#     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> # Where headers will be installed
# )

# Link with BLAS and LAPACK
target_link_libraries(slicot
    PUBLIC
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
)

# Install the library target
install(TARGETS slicot
    EXPORT SLICOTTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install C wrapper headers from include/ if enabled
if(SLICOT_BUILD_C_WRAPPERS)
    # Glob all headers from the main include directory
    file(GLOB C_WRAPPER_HEADERS "${CMAKE_SOURCE_DIR}/include/*.h")

    # Install the headers into include/slicot
    install(FILES ${C_WRAPPER_HEADERS}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/slicot
            COMPONENT Headers)
endif()
