# SLICOT Installation Guide

This guide describes how to build and install the SLICOT library using either the modern **CMake** build system (recommended) or the legacy **Makefile** system.

The essential source code and documentation for the **SLICOT** software is stored in the **SLICOT** Library root directory, **`slicot`**, and its subdirectories (**`benchmark_data`**, **`doc`**, **`examples`**, **`src`** and **`src_aux`**).

## Prerequisites

**Common:**

*   **C/C++ Compiler:** A compatible C/C++ compiler (needed by some dependencies or build tools).
*   **Fortran Compiler:** A Fortran compiler (e.g., `gfortran`, Intel Fortran `ifort`/`ifx`).
*   **BLAS and LAPACK:** Implementations of BLAS and LAPACK libraries. CMake will attempt to find these automatically. For Makefiles, you need to configure their paths in `make*.inc`. Common sources include:
    *   Intel MKL
    *   OpenBLAS
    *   System-provided libraries (e.g., via `apt`, `yum`, `brew`, MSYS2/MinGW packages)
    *   Netlib reference implementations (requires compilation).
    *   *Note:* For maximum efficiency, it is recommended to use machine-specific, optimized versions whenever possible.

**For CMake Build (Recommended):**

*   **CMake:** Version 3.15 or later.
*   **Build Tool:** A build tool supported by CMake (e.g., Ninja, Make, Visual Studio). Ninja is recommended and used in presets.

**For Legacy Makefile Build:**

*   **Make Tool:** A `make` utility (e.g., GNU Make on Linux/macOS, `nmake` on Windows with Visual Studio Command Prompt, or `mingw32-make` with MinGW).

## Building with CMake Presets (Recommended)

This project uses CMake Presets (`CMakePresets.json`) to simplify configuration and building across different platforms and compilers.

1.  **Choose a Preset:** Open the project in a CMake-aware IDE (like VS Code with CMake Tools, CLion) or use the command line. Ensure your chosen compiler toolchain (e.g., MinGW-w64 with gcc/g++/gfortran) is in your system's PATH. The available presets define common configurations:
    *   `windows-x64-debug` / `windows-x64-release` (Uses MinGW/Ninja by default)
    *   `linux-x64-debug` / `linux-x64-release` (Uses GCC/Ninja by default)
    *   `macos-x64-debug` / `macos-x64-release` (Uses Clang/gfortran/Ninja by default)

2.  **Configure:**
    *   **IDE:** Select the desired Configure Preset (e.g., `windows-x64-debug`). The IDE should run CMake automatically.
    *   **Command Line:**
        ```bash
        # Example for Windows Debug (ensure MinGW compilers are in PATH)
        cmake --preset windows-x64-debug

        # Example for Linux Debug
        # cmake --preset linux-x64-debug
        ```
        Replace `windows-x64-debug` with your chosen preset name. You can enable examples using `-DSLICOT_BUILD_EXAMPLES=ON`.

3.  **Build:**
    *   **IDE:** Select the corresponding Build Preset (e.g., `Build Windows x64 Debug (MinGW)`) or build the default target.
    *   **Command Line:**
        ```bash
        # Example for Windows Debug build
        cmake --build --preset windows-x64-debug-build

        # Example for Linux Debug build
        # cmake --build --preset linux-x64-debug-build
        ```
        Use the build preset name corresponding to your configure preset (usually adding `-build`).

4.  **Install (Optional but Recommended):** The build process creates library files in the build directory. To install them to the location specified by `CMAKE_INSTALL_PREFIX` in the preset (defaults to `./install/<preset-name>`), run the install step:
    *   **IDE:** Build the `INSTALL` target.
    *   **Command Line:**
        ```bash
        # Example for Windows Debug install using build preset
        cmake --build --preset windows-x64-debug-build --target install
        ```

5.  **Run Tests (if `SLICOT_BUILD_EXAMPLES=ON`):**
    *   **IDE:** Run tests using the IDE's test runner, selecting the appropriate Test Preset (e.g., `Test Windows x64 Debug (MinGW)`).
    *   **Command Line:**
        ```bash
        # Example for Windows Debug test using test preset
        ctest --preset windows-x64-debug-test
        ```

### CMake Customization (CMakeUserPresets.json)

If you need to specify a different compiler path, add custom compiler flags, or change other cache variables *without modifying the main presets file*, create a `CMakeUserPresets.json` file in the project root directory. This file allows you to inherit from the base presets in `CMakePresets.json` and override specific settings. See the provided `CMakeUserPresets.json` for examples.

### CMake Options

You can customize the build by setting CMake options during the configure step (e.g., using `-D<OPTION_NAME>=<VALUE>` on the command line or through IDE settings).

*   **`CMAKE_INSTALL_PREFIX`**: Specifies the directory where the library and documentation will be installed. Defaults are set in the presets (`install/<preset-name>`).
*   **`CMAKE_BUILD_TYPE`**: Specifies the build configuration (e.g., `Debug`, `Release`). Set by the presets.
*   **`CMAKE_C_COMPILER`**, **`CMAKE_CXX_COMPILER`**, **`CMAKE_Fortran_COMPILER`**: Specify the paths to the compilers. Defaults are set in the presets, but can be overridden (e.g., in `CMakeUserPresets.json`).
*   **`SLICOT_BUILD_EXAMPLES`**: Build the example programs and enable tests. Defaults to `OFF`. Set to `ON` to build and test the examples.
    ```bash
    # Example: Configure with examples enabled
    cmake --preset windows-x64-debug -DSLICOT_BUILD_EXAMPLES=ON
    ```
*   **`BUILD_SHARED_LIBS`**: Build SLICOT as a shared library (`.dll`, `.so`, `.dylib`) instead of a static library (`.lib`, `.a`). Defaults to `OFF`.

## Building with Legacy Makefiles

Template make files are provided to help building the **SLICOT** Library object file, and to link and run the available example programs calling the **SLICOT** Library routines.

1.  **Configure Makefiles:** In order to use these make files on a specific Unix-like or Windows platform, some changes might be needed in the files **`make*.inc`** and **`makefile*`** stored in the **SLICOT** (sub-)directories, **`slicot`**, **`examples`**, **`src`** and **`src_aux`**.
    *   The file named **`make.inc`** and the files **`makefile`** have been used on Windows platforms with Intel Fortran compilers.
    *   The files named **`make_Unix.inc`** and **`makefile_Unix`** are templates for Unix-like machines, including Linux, with gfortran compiler.
    *   Denote by <**slicotroot**> the path to the **`slicot`** directory (e.g., **`c:\slicot`**).
    *   The changes in **`make*.inc`** might define the specific machine (platform) identifier, the compiler, linker, and archiver flags, and the location and names of the **LAPACK** and **BLAS** libraries, which the program files should be linked to. Some details are given in the **`make*.inc`** files.

2.  **Build:** After performing the necessary changes, as suggested in the comments of the make files, the **SLICOT** library and examples can be generated automatically with the command:
    ```bash
    # On Unix-like systems (Linux, macOS) from the <slicotroot> directory
    make -f makefile_Unix

    # On Windows (e.g., using nmake from a Visual Studio Command Prompt)
    # from the <slicotroot> directory
    nmake -f makefile
    ```
    (Adjust `make`/`nmake` and the makefile name based on your system and chosen configuration).

3.  **Generated Files:** The first execution of **`(n)make`** will create the following files:
    *   The **SLICOT** Library object files **`*.o`** (Unix) or **`*.obj`** (Windows) in the subdirectories **`src`** and **`src_aux`**.
    *   The library files **`slicot.a`** and **`lpkaux.a`** (Unix), or **`slicot.lib`** and **`lpkaux.lib`** (Windows), in the directory **`slicot`**. The `lpkaux` libraries contain object files for LAPACK compatibility routines.
    *   The example programs object and executable files in the subdirectory **`examples`**.
    *   The files **`*.exa`**, with the results computed on the local machine, in the subdirectory **`examples`**.

4.  **Running Examples & Comparing Results:** The `(n)make all` target (default) in the `examples` directory typically runs the executables, redirecting input from `*.dat` files and output to `*.exa` files. These `*.exa` files may be compared with the reference results (`*.res`). Several types of differences could be noticed, including possible sign changes for some elements, or even different values in some columns and/or rows of the computed matrices (e.g., transformation matrices). This does not usually mean that the computed results are wrong.

5.  **Cleanup:** More details for executing other tasks, e.g., cleaning the object files (`make clean`) or object and executable files (`make cleanup`), are given in the `makefile*` files in the `slicot`, `src`, `src_aux`, and `examples` directories.

## Finding BLAS/LAPACK (General Notes)

Both build systems require BLAS and LAPACK.

*   **CMake:** Uses `find_package(BLAS)` and `find_package(LAPACK)`. Ensure they are installed and discoverable by CMake.
*   **Makefiles:** Require manual configuration of paths and library names in the `make*.inc` file.

Ensure the libraries are installed and discoverable.
*   On **Windows with MinGW**, you can install libraries like OpenBLAS via MSYS2/MinGW packages (e.g., `pacman -S mingw-w64-x86_64-openblas`). Alternatively, you can use a package manager like [vcpkg](https://vcpkg.io/) to install `openblas` or `mkl`. For CMake, integrate vcpkg by setting `-DCMAKE_TOOLCHAIN_FILE=[path/to/vcpkg]/scripts/buildsystems/vcpkg.cmake`. For Makefiles, you would manually set the library paths in `make.inc`. You can also download pre-built binaries and set environment variables like `OpenBLAS_HOME`.
*   On **Linux**, use your system package manager (e.g., `sudo apt install libopenblas-dev liblapack-dev` or `sudo yum install openblas-devel lapack-devel`).
*   On **macOS**, use Homebrew (e.g., `brew install openblas lapack`).

You might need to set environment variables (`MKLROOT`, `OpenBLAS_HOME`) or CMake variables (`-DBLAS_LIBRARIES=...`, `-DLAPACK_LIBRARIES=...`, `-DCMAKE_PREFIX_PATH=...`) if CMake cannot find the libraries automatically. For Makefiles, adjust the `BLASLIB` and `LAPACKLIB` variables in `make*.inc`.

## Using the Installed Library (CMake Build)

After installation using the CMake build system, you can use the SLICOT library in other CMake projects:

```cmake
# Find the installed SLICOT package
# Ensure CMAKE_PREFIX_PATH includes the SLICOT install directory if needed
find_package(SLICOT REQUIRED)

# Link your target against the imported SLICOT library target
target_link_libraries(your_target PRIVATE SLICOT::slicot)
```

The legacy Makefile build does not provide a standard installation or CMake package configuration. You would typically need to manually specify the path to `slicot.lib`/`slicot.a` and link against it and its dependencies (BLAS, LAPACK) in your own project's build system.

