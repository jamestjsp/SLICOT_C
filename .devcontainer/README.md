# SLICOT C Wrapper Dev Container

This directory contains the development container configuration for the SLICOT C Wrapper project. The dev container provides a consistent, reproducible development environment with all necessary tools and dependencies pre-installed.

## What's Included

### Compilers and Build Tools
- **GCC/G++ 13.x** - C/C++ compilers
- **GFortran 13.x** - Fortran 77/90/95 compiler
- **CMake 3.22+** - Modern build system
- **Ninja** - Fast parallel build system
- **Make** - Traditional build system

### Libraries
- **OpenBLAS** - Optimized BLAS implementation
- **LAPACK** - Linear algebra package
- **GoogleTest** - C++ testing framework (fetched by CMake)

### Development Tools
- **Git & Git LFS** - Version control
- **GDB** - GNU debugger
- **Valgrind** - Memory debugging and profiling
- **Python 3** - For scripts and utilities (includes NumPy, SciPy, Matplotlib, Pandas)

### VS Code Extensions
- C/C++ Extension Pack (IntelliSense, debugging, etc.)
- CMake Tools (CMake integration)
- Modern Fortran (Fortran language support)
- GitLens (enhanced Git integration)
- Code Spell Checker
- Error Lens (inline error display)
- Todo Tree (TODO/FIXME highlighting)

## Quick Start

### 1. Prerequisites
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop)
- Install [Visual Studio Code](https://code.visualstudio.com/)
- Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

### 2. Open in Container
1. Open this project folder in VS Code
2. Press `F1` and select **"Dev Containers: Reopen in Container"**
3. Wait for the container to build (first time takes 5-10 minutes)
4. The `postCreate.sh` script will automatically configure CMake with the `linux-x64-debug` preset

### 3. Build and Test
Once the container is ready, you can build and test:

```bash
# Build the project (Debug configuration)
cmake --build --preset linux-x64-debug-build

# Run all tests
ctest --preset linux-x64-debug-test

# Run a specific test
./build/linux-x64-debug/tests/slicot_c_wrapper_tests --gtest_filter="AB01MD*"
```

## Using Different Build Presets

The project supports multiple CMake presets. In the dev container, use the Linux presets:

### Debug Build (Default)
```bash
cmake --preset linux-x64-debug
cmake --build --preset linux-x64-debug-build
ctest --preset linux-x64-debug-test
```

### Release Build
```bash
cmake --preset linux-x64-release
cmake --build --preset linux-x64-release-build
ctest --preset linux-x64-release-test
```

### Build with Examples
```bash
cmake --preset linux-x64-debug -DSLICOT_BUILD_EXAMPLES=ON
cmake --build --preset linux-x64-debug-build
```

### Build Options
You can customize the build by setting CMake cache variables:

```bash
cmake --preset linux-x64-debug \
  -DSLICOT_BUILD_TESTS=ON \
  -DSLICOT_BUILD_C_WRAPPERS=ON \
  -DSLICOT_BUILD_EXAMPLES=ON \
  -DSLICOT_BUILD_SHARED_LIBS=ON
```

Available options:
- `SLICOT_BUILD_TESTS` - Build GoogleTest test suite (default: ON)
- `SLICOT_BUILD_C_WRAPPERS` - Build C wrapper layer (default: ON)
- `SLICOT_BUILD_EXAMPLES` - Build example programs (default: OFF)
- `SLICOT_BUILD_SHARED_LIBS` - Build shared libraries (default: ON)
- `SLICOT_USE_ILP64` - Use 64-bit integers (default: OFF)

## Directory Structure

The container creates volumes for build artifacts to improve performance:

```
/workspaces/SLICOT_C/          # Your project root
├── build/                      # Build output (Docker volume)
│   └── linux-x64-debug/        # Debug build artifacts
├── install/                    # Installation output (Docker volume)
│   └── linux-x64-debug/        # Debug install artifacts
└── ...                         # Source files
```

## Running Tests

### All Tests
```bash
ctest --preset linux-x64-debug-test
```

### Specific Test Suite
```bash
./build/linux-x64-debug/tests/slicot_c_wrapper_tests --gtest_filter="AB01MD*"
```

### Verbose Test Output
```bash
ctest --preset linux-x64-debug-test --verbose
```

### Failed Tests Only
```bash
ctest --preset linux-x64-debug-test --rerun-failed --output-on-failure
```

## Debugging

### Using GDB
```bash
# Build with debug symbols (already enabled in debug preset)
cmake --build --preset linux-x64-debug-build

# Run with GDB
gdb ./build/linux-x64-debug/tests/slicot_c_wrapper_tests

# In GDB:
(gdb) break main
(gdb) run --gtest_filter="AB01MD*"
```

### Using VS Code Debugger
1. Set breakpoints in your code
2. Press `F5` or use the Debug panel
3. Select the appropriate launch configuration (pre-configured in `.vscode/launch.json` if available)

### Memory Debugging with Valgrind
```bash
valgrind --leak-check=full \
  ./build/linux-x64-debug/tests/slicot_c_wrapper_tests \
  --gtest_filter="AB01MD*"
```

## Installing the Library

To install SLICOT to the system (within the container):

```bash
cmake --build --preset linux-x64-debug-build --target install
```

Default install location: `install/linux-x64-debug/`

To change the install prefix:
```bash
cmake --preset linux-x64-debug \
  -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build --preset linux-x64-debug-build --target install
```

## Customizing the Container

### Modifying the Dockerfile
1. Edit `.devcontainer/Dockerfile`
2. Rebuild the container:
   - Press `F1`
   - Select **"Dev Containers: Rebuild Container"**

### Modifying VS Code Settings
1. Edit `.devcontainer/devcontainer.json`
2. Reload the window:
   - Press `F1`
   - Select **"Developer: Reload Window"**

### Adding Extensions
Add extension IDs to the `extensions` array in `devcontainer.json`:

```json
"customizations": {
  "vscode": {
    "extensions": [
      "ms-vscode.cpptools",
      "your-extension-id-here"
    ]
  }
}
```

## Performance Optimization

The dev container uses Docker volumes for `build/` and `install/` directories to improve I/O performance. If you experience slow builds:

1. **Use Ninja instead of Make** (already configured)
2. **Parallel builds**: `cmake --build --preset linux-x64-debug-build -j$(nproc)`
3. **Incremental builds**: Only rebuild what changed

## Troubleshooting

### Container won't start
- Ensure Docker Desktop is running
- Check Docker Desktop has sufficient resources (4GB+ RAM recommended)
- Try rebuilding: `F1` → "Dev Containers: Rebuild Container"

### CMake configuration fails
```bash
# Clean build directory
rm -rf build/linux-x64-debug

# Reconfigure
cmake --preset linux-x64-debug
```

### Tests fail
- Check `tests/Known_Issues.txt` for known failures
- Ensure you're using the debug preset for testing
- Run with verbose output: `ctest --preset linux-x64-debug-test --verbose`

### Git issues
If you see "detected dubious ownership" errors:
```bash
git config --global --add safe.directory /workspaces/SLICOT_C
```

## Compiler Versions

The container uses GCC/G++/GFortran 13.x by default. To verify:

```bash
gcc --version
g++ --version
gfortran --version
```

## BLAS/LAPACK Configuration

The container uses OpenBLAS and LAPACK from Ubuntu repositories. CMake should automatically detect them. To verify:

```bash
dpkg -l | grep -E 'openblas|lapack'
```

## Python Environment

Python 3 is installed with the following packages:
- NumPy
- SciPy
- Matplotlib
- Pandas

To install additional packages:
```bash
pip3 install --user package-name
```

## Additional Resources

- **Project Documentation**: See `CLAUDE.md` for detailed project guidelines
- **CMake Presets**: See `CMakePresets.json` for all available presets
- **Build Instructions**: See `README.md` and `Installation.md`
- **Test Data**: See `tests/Readme.md` for test data format guidelines

## Container Architecture

- **Base Image**: Ubuntu 24.04 LTS
- **User**: `vscode` (non-root, UID 1000)
- **Workspace**: `/workspaces/SLICOT_C`
- **Build Dir**: `/workspaces/SLICOT_C/build` (Docker volume)
- **Install Dir**: `/workspaces/SLICOT_C/install` (Docker volume)

## License

This dev container configuration is part of the SLICOT C wrapper project. See the main project LICENSE file for details.
