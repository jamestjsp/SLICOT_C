# Build SLICOT examples and tests
include(GNUInstallDirs)
enable_testing()

# Find dependencies (already found in parent scope, but good practice)
# SLICOT library target 'slicot' is available from parent scope

# Find all example source files (convention: T*.f)
file(GLOB EXAMPLE_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "T*.f")
# Find all data files (convention: *.dat)
file(GLOB DATA_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.dat")

# Process each example source file
foreach(EXAMPLE_SRC ${EXAMPLE_SOURCES})
    # Get the base name (e.g., TAB01MD from TAB01MD.f)
    string(REGEX REPLACE "\\.f$" "" EXAMPLE_NAME ${EXAMPLE_SRC})

    # Add the executable target
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SRC})
    target_link_libraries(${EXAMPLE_NAME} PRIVATE slicot ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    set_target_properties(${EXAMPLE_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

    # Check if a corresponding data file exists
    # T<name>.f uses <name>.dat
    # TT<name>.f uses T<name>.dat
    if(EXAMPLE_NAME MATCHES "^TT")
      # Example: TTB01ID.f -> TB01ID.dat
      string(REGEX REPLACE "^TT" "T" DATA_FILE_BASE ${EXAMPLE_NAME})
    else()
      # Example: TAB01MD.f -> AB01MD.dat
      string(REGEX REPLACE "^T" "" DATA_FILE_BASE ${EXAMPLE_NAME})
    endif()

    set(DATA_FILE "${DATA_FILE_BASE}.dat")
    set(RESULT_FILE "${DATA_FILE_BASE}.exa") # Expected output file name matches data file base
    set(GOLD_RESULT_FILE "${DATA_FILE_BASE}.res") # Reference result file name matches data file base
    set(TEST_NAME "${DATA_FILE_BASE}_test") # Test name matches data file base

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${DATA_FILE}")
        # Add a test for this example
        add_test(
            NAME ${TEST_NAME}
            COMMAND $<TARGET_FILE:${EXAMPLE_NAME}>
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} # Run test in build dir
        )

        # Configure the test to use the data file as input and check output
        set_tests_properties(${TEST_NAME} PROPERTIES
            ATTACHED_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${DATA_FILE}" # Make data file available
            # Redirect input from the data file
            COMMAND_MODIFICATION "< \"${CMAKE_CURRENT_SOURCE_DIR}/${DATA_FILE}\" > \"${RESULT_FILE}\""
        )

        # Optionally, compare the output with a reference result if it exists
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${GOLD_RESULT_FILE}")
             set_tests_properties(${TEST_NAME} PROPERTIES
                 PASS_REGULAR_EXPRESSION "Test passed" # Default, adjust if needed
                 # Add command to compare output file with reference file after execution
                 # Example using cmake -E compare_files (basic diff):
                 # POST_COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol "${CMAKE_CURRENT_BINARY_DIR}/${RESULT_FILE}" "${CMAKE_CURRENT_SOURCE_DIR}/${GOLD_RESULT_FILE}"
             )
             # Add the reference file as an attached file for reference
             # Need to collect attached files properly if both data and gold exist
             set(test_attached_files "${CMAKE_CURRENT_SOURCE_DIR}/${DATA_FILE}") # Start with data file
             list(APPEND test_attached_files "${CMAKE_CURRENT_SOURCE_DIR}/${GOLD_RESULT_FILE}")
             set_property(TEST ${TEST_NAME} PROPERTY ATTACHED_FILES ${test_attached_files})
        endif()

    else()
        message(STATUS "Data file ${DATA_FILE} not found for example ${EXAMPLE_NAME}, test not created.")
    endif()

    # Install the example executable (optional)
    # install(TARGETS ${EXAMPLE_NAME}
    #         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/examples COMPONENT Examples)

endforeach()

# Install example source files, data files, and result files
file(GLOB EXAMPLE_F_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.f")
file(GLOB EXAMPLE_DAT_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.dat")
file(GLOB EXAMPLE_RES_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.res")

# Check if file lists are not empty before installing
if(EXAMPLE_F_FILES)
    install(FILES ${EXAMPLE_F_FILES}
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
            COMPONENT ExamplesDoc)
endif()
if(EXAMPLE_DAT_FILES)
    install(FILES ${EXAMPLE_DAT_FILES}
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
            COMPONENT ExamplesDoc)
endif()
if(EXAMPLE_RES_FILES)
    install(FILES ${EXAMPLE_RES_FILES}
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
            COMPONENT ExamplesDoc)
endif()
