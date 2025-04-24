# CMakeLists.txt for SLICOT C Wrappers

include(${CMAKE_SOURCE_DIR}/cmake/slicot_sources.cmake) # Include the source list definition

# Create an object library for the C wrappers
add_library(slicot_c_wrapper_obj OBJECT ${SLICOT_C_WRAPPER_SOURCES}) # Use the defined list

# Make the object library part of the export set
set_property(TARGET slicot_c_wrapper_obj PROPERTY EXPORT_NAME SLICOT::slicot_c_wrapper_obj)
set_property(TARGET slicot_c_wrapper_obj PROPERTY POSITION_INDEPENDENT_CODE ON)

# Define macros for export/import or static build based on SLICOT_BUILD_SHARED_LIBS
# Also set WINDOWS_EXPORT_ALL_SYMBOLS for shared Windows builds
if(SLICOT_BUILD_SHARED_LIBS)
    # Define SLICOT_C_WRAPPER_EXPORTS when building the shared library
    target_compile_definitions(slicot_c_wrapper_obj PRIVATE SLICOT_C_WRAPPER_EXPORTS)
    # Automatically export symbols from object files for Windows shared builds
    if(WIN32)
        set_target_properties(slicot_c_wrapper_obj PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
        message(STATUS "Setting WINDOWS_EXPORT_ALL_SYMBOLS for slicot_c_wrapper_obj.")
    endif()
else()
    # Define SLICOT_STATIC when building a static library
    target_compile_definitions(slicot_c_wrapper_obj PUBLIC SLICOT_STATIC)
endif()

# Set target properties if needed, e.g., include directories
target_include_directories(slicot_c_wrapper_obj PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include # Include directory for building
)

# Add dependencies if necessary, e.g., if wrappers depend on slicot headers
# target_link_libraries(slicot_c_wrapper_obj INTERFACE slicot) # Example if needed

message(STATUS "Configuring SLICOT C wrapper object library.")
