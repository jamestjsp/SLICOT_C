# Building SLICOT on Windows

## 🚀 Quick Start for Junior Developers

**New to development on Windows? Start here!** This simplified guide will get you up and running quickly.

### Step 1: Install Compiler and Toolchain

1. **Download w64devkit** (all-in-one GCC toolchain):
   - Go to https://github.com/skeeto/w64devkit/releases
   - Download the latest release (it's a self-extracting zip)
   - Extract to `C:\w64devkit` (or your preferred location)

2. **Add to PATH**:
   - Press `Windows + R`, type `sysdm.cpl`, and press Enter
   - Go to **Advanced** tab → **Environment Variables**
   - Under "User variables" or "System variables", select `Path` → Click **Edit**
   - Click **New** and add: `C:\w64devkit\bin` (or wherever you extracted it)
   - Click **OK** on all windows

3. **Set compiler environment variables**:
   - In the same Environment Variables window, click **New** for each:
     - Variable name: `CC`, Value: `gcc`
     - Variable name: `FC`, Value: `gfortran`
     - Variable name: `CXX`, Value: `g++`
   - Click **OK** to save

4. **Verify installation**:
   - Open a **new** PowerShell window (important: must be new!)
   - Type these commands to check:
     ```powershell
     gcc --version
     gfortran --version
     g++ --version
     ```
   - You should see version information for each

### Step 2: Install CMake and Ninja

Open PowerShell **as Administrator** and run:

```powershell
winget install Kitware.CMake --version 3.30.1
winget install Ninja-build.Ninja
```

Close and reopen PowerShell, then verify:
```powershell
cmake --version
ninja --version
```

### Step 3: Install vcpkg (Package Manager)

vcpkg manages libraries like BLAS and LAPACK that SLICOT needs.

```powershell
# Clone vcpkg
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
cd C:\vcpkg

# Bootstrap (install) vcpkg
.\bootstrap-vcpkg.bat

```

**Set VCPKG_ROOT environment variable**:
- Go back to Environment Variables (Step 1.2)
- Click **New** under User variables
- Variable name: `VCPKG_ROOT`, Value: `C:\vcpkg`
- Also add `C:\vcpkg` to your PATH
- Click **OK**

### Step 4: Install Visual Studio Code

1. **Download VS Code**: https://code.visualstudio.com/
2. **Install it** (default options are fine)
3. **Install extensions**:
   - Open VS Code
   - Press `Ctrl+Shift+X` (opens Extensions panel)
   - Search for and install:
     - **C/C++** (by Microsoft) - ID: `ms-vscode.cpptools`
     - **CMake Tools** (by Microsoft) - ID: `ms-vscode.cmake-tools`

### Step 5: Build SLICOT (Using VS Code GUI)

1. **Open the project**:
   - In VS Code: File → Open Folder
   - Select the `SLICOT_C` root folder
   - Wait a moment for the CMake Tools extension to activate

2. **Select a Configure Preset** (bottom status bar):
   - Look at the bottom status bar in VS Code
   - Click on **"No Configure Preset Selected"** (or the current preset name)
   - From the dropdown, choose: `windows-x64-release-mingw`
   - The project will automatically configure (you'll see progress in the Output panel)

3. **Select a Build Preset** (bottom status bar):
   - Click on **"No Build Preset Selected"** (or the current preset name)
   - Choose: `windows-x64-release-mingw`

4. **Build the project**:
   - Click the **Build** button (hammer icon) in the bottom status bar
   - OR press `F7`
   - Watch the build progress in the Output panel
   - Wait for "Build finished" message

5. **Run tests** (optional):
   - Click the **Run CTest** button (play icon with checkmark) in the bottom status bar
   - View test results in the Output panel

**Visual Guide - Bottom Status Bar:**
```
[CMake Preset] [Build Preset] [🔨 Build] [▶ Run] [🐛 Debug] [✓ Run CTest]
```

**Alternative: Command Palette Method**
- Press `Ctrl+Shift+P` and type:
  - `CMake: Select Configure Preset` → Choose `windows-x64-release-mingw`
  - `CMake: Configure`
  - `CMake: Build`

**Alternative: Command Line Build**
```powershell
# Navigate to project folder
cd C:\Users\JOSEPHJ\Documents\dev\SLICOT_C

# Configure
cmake --preset windows-x64-release-mingw

# Build
cmake --build build/windows-x64-release-mingw
```

### Common Issues

**"Command not found" errors:**
- Make sure you restarted PowerShell after adding things to PATH
- Check Environment Variables are set correctly

**CMake can't find compiler:**
- Verify environment variables: `echo $env:CC` should show `gcc`
- Restart VS Code after setting environment variables

**vcpkg errors:**
- Make sure `VCPKG_ROOT` is set and you restarted PowerShell/VS Code

### Need More Details?

See the full guide below for advanced options and troubleshooting.

---

## 📚 Full Build Guide

This guide provides step-by-step instructions for building SLICOT on Windows using different toolchains and vcpkg for dependency management.

## Prerequisites

- **CMake** (version 3.30.1 recommended) - [Download](https://cmake.org/download/) or install via `winget install Kitware.CMake --version 3.30.1`
- **Git** - [Download](https://git-scm.com/downloads)
- **Ninja** build system (recommended) - [Download](https://github.com/ninja-build/ninja/releases) or install via `winget install Ninja-build.Ninja`
- A C/C++ compiler and Fortran compiler (see toolchain options below)

## Toolchain Options

Choose one of the following toolchains:

### Option 1: Intel oneAPI (Recommended for Intel MKL)
- Intel oneAPI Base Toolkit (includes Intel MKL)
- Intel oneAPI HPC Toolkit (includes Intel Fortran Compiler)
- [Download Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)

### Option 2: MinGW-w64 with vcpkg (Recommended for open-source stack)
- MinGW-w64 GCC/Gfortran - [Download](https://www.mingw-w64.org/)
- Use vcpkg to install BLAS/LAPACK (OpenBLAS or LAPACK Reference)

### Option 3: MSVC + Intel Fortran
- Visual Studio 2019 or later with C++ support
- Intel Fortran Compiler (ifx/ifort)

## Setting Up vcpkg

vcpkg is Microsoft's package manager for C/C++ libraries. It simplifies dependency management for BLAS and LAPACK.

### 1. Clone and Bootstrap vcpkg

Open PowerShell and run:

```powershell
# Clone vcpkg repository
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
cd C:\vcpkg

# Bootstrap vcpkg
.\bootstrap-vcpkg.bat
```

### 2. Set Environment Variables

Set the `VCPKG_ROOT` environment variable (required for CMake integration):

**Temporary (current session only):**
```powershell
$env:VCPKG_ROOT = "C:\vcpkg"
$env:PATH = "$env:VCPKG_ROOT;$env:PATH"
```

**Permanent (recommended):**
```powershell
# Set system environment variable (requires admin privileges)
[System.Environment]::SetEnvironmentVariable('VCPKG_ROOT', 'C:\vcpkg', 'User')
[System.Environment]::SetEnvironmentVariable('PATH', "C:\vcpkg;$env:PATH", 'User')
```

Or manually:
1. Open Windows Settings → System → About → Advanced system settings
2. Click "Environment Variables"
3. Add a new User variable: `VCPKG_ROOT` = `C:\vcpkg`
4. Add `C:\vcpkg` to your PATH

### 3. Install BLAS/LAPACK Dependencies

Choose one of the following BLAS/LAPACK implementations:

**OpenBLAS (recommended for MinGW):**
```powershell
vcpkg install openblas:x64-mingw-dynamic
vcpkg install lapack:x64-mingw-dynamic
```

**LAPACK Reference Implementation:**
```powershell
vcpkg install lapack-reference:x64-mingw-dynamic
```

**For other triplets:**
- `x64-mingw-static` - Static libraries for MinGW
- `x64-windows` - Dynamic libraries for MSVC
- `x64-windows-static` - Static libraries for MSVC

### 4. Verify Installation

```powershell
vcpkg list
```

You should see your installed packages listed.

## Building SLICOT

### Step 1: Clone the SLICOT Repository

```powershell
git clone https://github.com/jamestjsp/SLICOT_C.git
cd SLICOT_C
```

### Step 2: Configure with CMake

SLICOT uses CMake presets for easy configuration. Choose a preset based on your toolchain:

#### Using MinGW with vcpkg (Debug):

```powershell
cmake --preset=windows-x64-debug-mingw
```

#### Using MinGW with vcpkg (Release):

```powershell
cmake --preset=windows-x64-release-mingw
```

#### Using Intel Compilers (Debug):

First, initialize the Intel compiler environment:
```powershell
# Adjust path based on your Intel oneAPI installation
& "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```

Then configure:
```powershell
cmake --preset=windows-x64-debug-intel
```

#### Using Intel Compilers (Release):

```powershell
& "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cmake --preset=windows-x64-release-intel
```

### Step 3: Build the Project

After configuration, build the project:

```powershell
cmake --build build/<preset-name>
```

For example:
```powershell
# For MinGW Debug
cmake --build build/windows-x64-debug-mingw

# For MinGW Release
cmake --build build/windows-x64-release-mingw

# For Intel Debug
cmake --build build/windows-x64-debug-intel
```

#### Build with Parallel Jobs (faster):

```powershell
cmake --build build/<preset-name> --parallel 8
```

### Step 4: Run Tests (Optional)

```powershell
cd build/<preset-name>
ctest --output-on-failure
```

### Step 5: Install (Optional)

```powershell
cmake --install build/<preset-name>
```

This installs SLICOT to the directory specified in the preset (typically `install/<preset-name>`).

## CMake Build Options

You can customize the build with these options (set in CMakeLists.txt or via `-D` flag):

- `SLICOT_BUILD_SHARED_LIBS` - Build shared libraries (ON by default)
- `SLICOT_BUILD_TESTS` - Build C wrapper tests with GoogleTest (ON by default)
- `SLICOT_BUILD_EXAMPLES` - Build example programs (OFF by default)
- `SLICOT_BUILD_C_WRAPPERS` - Build C wrapper library (ON by default)
- `SLICOT_USE_VCPKG` - Use vcpkg for dependencies (ON for MinGW presets)
- `SLICOT_LINK_AUX_WRAPPERS` - Link auxiliary LAPACK wrappers (ON by default)

Example with custom options:
```powershell
cmake --preset=windows-x64-release-mingw -DSLICOT_BUILD_EXAMPLES=ON
```

## Troubleshooting

### vcpkg not found

**Error:** `CMAKE_TOOLCHAIN_FILE not found` or `vcpkg.cmake not found`

**Solution:** Ensure `VCPKG_ROOT` environment variable is set correctly and restart your terminal/IDE.

### BLAS/LAPACK not found

**Error:** `Could NOT find BLAS` or `Could NOT find LAPACK`

**Solution:**
1. Verify packages are installed: `vcpkg list`
2. Ensure you're using a preset with `SLICOT_USE_VCPKG=ON` (MinGW presets)
3. Check that the triplet matches your compiler (e.g., `x64-mingw-dynamic` for MinGW)

### Compiler not found

**Error:** `CMAKE_Fortran_COMPILER not found`

**Solution:**
- For Intel: Run `setvars.bat` before CMake configuration
- For MinGW: Ensure MinGW bin directory is in your PATH
- Verify compiler: `gfortran --version` or `ifx --version`

### Ninja not found

**Error:** `CMake Error: CMake was unable to find a build program corresponding to "Ninja"`

**Solution:**
1. Download Ninja from https://github.com/ninja-build/ninja/releases
2. Add Ninja to your PATH, or
3. Change generator in preset to `"Visual Studio 17 2022"` or `"Unix Makefiles"`

### Mixed string length arguments (Intel Fortran)

**Error:** String length argument mismatch

**Solution:** SLICOT CMakeLists.txt already sets the correct flags (`/iface:nomixed_str_len_arg`). Ensure you're using the Intel presets.

## Additional Resources

- [vcpkg Documentation](https://learn.microsoft.com/en-us/vcpkg/get_started/get-started)
- [CMake Documentation](https://cmake.org/documentation/)
- [SLICOT Documentation](http://slicot.org/)
- [Intel oneAPI Documentation](https://www.intel.com/content/www/us/en/developer/tools/oneapi/documentation.html)

## Quick Start Example

Here's a complete example for MinGW with vcpkg:

```powershell
# 1. Set up vcpkg
git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
cd C:\vcpkg
.\bootstrap-vcpkg.bat
$env:VCPKG_ROOT = "C:\vcpkg"

# 2. Install dependencies
vcpkg install openblas:x64-mingw-dynamic lapack:x64-mingw-dynamic

# 3. Clone and build SLICOT
cd C:\Users\YourUsername\Documents
git clone https://github.com/jamestjsp/SLICOT_C.git
cd SLICOT_C

# 4. Configure and build
cmake --preset=windows-x64-release-mingw
cmake --build build/windows-x64-release-mingw --parallel

# 5. Run tests
cd build/windows-x64-release-mingw
ctest --output-on-failure
```

## Notes

- The MinGW presets automatically use vcpkg via `CMAKE_TOOLCHAIN_FILE`
- Intel presets use Intel MKL directly (no vcpkg needed)
- CMake presets are defined in `CMakePresets.json`
- User-specific settings can be added to `CMakeUserPresets.json` (not tracked by git)

---

For more detailed information, see `Installation.md` and `README.md` in the project root.
