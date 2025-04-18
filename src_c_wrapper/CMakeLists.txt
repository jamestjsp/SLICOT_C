# CMakeLists.txt for SLICOT C Wrappers

# Find all C source files in the current directory
file(GLOB C_WRAPPER_SOURCES "*.c")

# Create an object library for the C wrappers
add_library(slicot_c_wrapper_obj OBJECT ${C_WRAPPER_SOURCES})

# Define macros for export/import or static build based on SLICOT_BUILD_SHARED_LIBS
if(SLICOT_BUILD_SHARED_LIBS)
    # Define SLICOT_C_WRAPPER_EXPORTS when building the shared library
    target_compile_definitions(slicot_c_wrapper_obj PRIVATE SLICOT_C_WRAPPER_EXPORTS)
else()
    # Define SLICOT_STATIC when building a static library
    target_compile_definitions(slicot_c_wrapper_obj PRIVATE SLICOT_STATIC)
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
