#!/bin/bash
# Post-creation script for SLICOT C Wrapper dev container
# This script runs after the container is created to initialize the development environment

set -e  # Exit on error

echo "=========================================="
echo "SLICOT C Wrapper Development Container"
echo "Post-Create Initialization"
echo "=========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the workspace directory
WORKSPACE_DIR="${CONTAINERWORKSPACEFOLDER:-/workspaces/SLICOT_C}"
cd "$WORKSPACE_DIR"

echo -e "${BLUE}Current directory: $(pwd)${NC}"
echo ""

# Configure git safe directory
echo -e "${YELLOW}Configuring git safe directory...${NC}"
git config --global --add safe.directory "$WORKSPACE_DIR"
echo -e "${GREEN}Git safe directory configured.${NC}"
echo ""

# Ensure uv CLI is available
echo -e "${YELLOW}Ensuring uv CLI is installed...${NC}"
if ! command -v uv >/dev/null 2>&1; then
    if curl -LsSf https://astral.sh/uv/install.sh | sh -s -- --yes; then
        echo -e "${GREEN}uv installed successfully.${NC}"
    else
        echo -e "${YELLOW}uv installation failed.${NC}"
    fi
else
    echo -e "${GREEN}uv already installed.${NC}"
fi

# Make sure the local bin directory is on PATH for the current session and future shells
if [ -d "$HOME/.local/bin" ]; then
    export PATH="$HOME/.local/bin:$PATH"
    touch "$HOME/.bashrc"
    if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc"; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
    fi
fi
echo ""

# Initialize and sync uv project environment
echo -e "${YELLOW}Syncing uv project environment...${NC}"
if command -v uv >/dev/null 2>&1; then
    if [ ! -f "pyproject.toml" ]; then
        echo -e "${BLUE}pyproject.toml not found; initializing minimal uv project...${NC}"
        uv init --bare --no-readme --no-pin-python --name slicot-c-tooling . || true
    fi

    if uv sync --frozen; then
        echo -e "${GREEN}uv environment synchronized.${NC}"
    else
        echo -e "${YELLOW}uv sync failed; continuing without project environment.${NC}"
    fi
else
    echo -e "${YELLOW}Skipping uv sync because uv is unavailable.${NC}"
fi
echo ""

# Ensure the project virtual environment binaries are on PATH
if [ -d ".venv/bin" ]; then
    VENV_BIN_DIR="$(pwd)/.venv/bin"
    export PATH="$VENV_BIN_DIR:$PATH"
    touch "$HOME/.bashrc"
    if ! grep -q "export PATH=\"$VENV_BIN_DIR:\$PATH\"" "$HOME/.bashrc"; then
        echo "export PATH=\"$VENV_BIN_DIR:\$PATH\"" >> "$HOME/.bashrc"
    fi
fi
echo ""

# Verify fortls availability via uv-managed environment
echo -e "${YELLOW}Verifying fortls language server...${NC}"
if command -v uv >/dev/null 2>&1; then
    if uv run -- fortls --version; then
        echo -e "${GREEN}fortls is ready for use.${NC}"
    else
        echo -e "${YELLOW}fortls verification failed; ensure the uv environment is configured correctly.${NC}"
    fi
else
    echo -e "${YELLOW}Skipping fortls verification because uv is unavailable.${NC}"
fi
echo ""

# Install Claude Code CLI
echo -e "${YELLOW}Installing Claude Code CLI...${NC}"
if npm install -g @anthropic-ai/claude-code 2>/dev/null; then
    echo -e "${GREEN}Claude Code CLI installed successfully!${NC}"
    claude --version
else
    echo -e "${YELLOW}Claude Code CLI installation skipped or failed.${NC}"
fi
echo ""

# Display compiler versions
echo -e "${YELLOW}Installed compiler versions:${NC}"
echo "-----------------------------------"
gcc --version | head -n1
g++ --version | head -n1
gfortran --version | head -n1
cmake --version | head -n1
ninja --version
echo ""

# Display BLAS/LAPACK info
echo -e "${YELLOW}BLAS/LAPACK libraries:${NC}"
echo "-----------------------------------"
dpkg -l | grep -E 'openblas|lapack|blas' | awk '{print $2, $3}'
echo ""

# Run initial CMake configure with linux-x64-debug preset
echo -e "${YELLOW}Running initial CMake configuration...${NC}"
echo "Preset: linux-x64-debug"
echo "-----------------------------------"

if cmake --preset linux-x64-debug; then
    echo -e "${GREEN}CMake configuration successful!${NC}"
    echo ""

    # Optionally build the project (commented out by default to save time)
    # Uncomment the following lines to build automatically on container creation
    # echo -e "${YELLOW}Building project...${NC}"
    # if cmake --build --preset linux-x64-debug-build -j$(nproc); then
    #     echo -e "${GREEN}Build successful!${NC}"
    # else
    #     echo -e "${YELLOW}Build failed. You can build manually later.${NC}"
    # fi
else
    echo -e "${YELLOW}CMake configuration failed. You can configure manually later.${NC}"
fi

echo ""
echo "=========================================="
echo -e "${GREEN}Container initialization complete!${NC}"
echo "=========================================="
echo ""
echo -e "${BLUE}Quick Start Guide:${NC}"
echo "-----------------------------------"
echo ""
echo "1. Build the project (Debug):"
echo "   $ cmake --build --preset linux-x64-debug-build"
echo ""
echo "2. Run tests:"
echo "   $ ctest --preset linux-x64-debug-test"
echo ""
echo "3. Run specific test:"
echo "   $ ./build/linux-x64-debug/tests/slicot_c_wrapper_tests --gtest_filter=\"AB01MD*\""
echo ""
echo "4. Switch to Release build:"
echo "   $ cmake --preset linux-x64-release"
echo "   $ cmake --build --preset linux-x64-release-build"
echo ""
echo "5. Build with examples:"
echo "   $ cmake --preset linux-x64-debug -DSLICOT_BUILD_EXAMPLES=ON"
echo "   $ cmake --build --preset linux-x64-debug-build"
echo ""
echo "6. Install to system:"
echo "   $ cmake --build --preset linux-x64-debug-build --target install"
echo ""
echo "7. Use Claude Code CLI:"
echo "   $ claude"
echo ""
echo -e "${YELLOW}Available CMake Presets:${NC}"
echo "  - linux-x64-debug     (Default, configured)"
echo "  - linux-x64-release"
echo ""
echo -e "${YELLOW}Build Directory:${NC} build/linux-x64-debug"
echo -e "${YELLOW}Install Directory:${NC} install/linux-x64-debug"
echo ""
echo -e "${BLUE}For more information, see:${NC}"
echo "  - CLAUDE.md (AI assistant guide)"
echo "  - README.md (project overview)"
echo "  - .devcontainer/README.md (container-specific docs)"
echo ""
echo "Happy coding!"
echo ""
