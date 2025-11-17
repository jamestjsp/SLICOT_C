# SLICOT C Wrapper Development Container
# Based on Ubuntu 24.04 LTS with complete build toolchain for Fortran/C/C++ development

FROM ubuntu:24.04

# Avoid interactive prompts during package installation
ARG DEBIAN_FRONTEND=noninteractive

# Set timezone to avoid tzdata prompts
ENV TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install base system dependencies and build tools
RUN apt-get update && apt-get install -y \
    # Essential build tools
    build-essential \
    # LLVM/Clang/Flang compilers (version 20 - latest)
    llvm-20 \
    llvm-20-dev \
    llvm-20-tools \
    clang-20 \
    flang-20 \
    libflang-20-dev \
    libc++-20-dev \
    libc++abi-20-dev \
    # Build system tools
    cmake \
    ninja-build \
    make \
    pkg-config \
    # Version control
    git \
    git-lfs \
    # BLAS/LAPACK libraries
    libopenblas-dev \
    liblapack-dev \
    libblas-dev \
    # Utilities
    wget \
    curl \
    zip \
    unzip \
    tar \
    gzip \
    ca-certificates \
    gnupg \
    software-properties-common \
    # Debugging tools
    gdb \
    lldb-20 \
    valgrind \
    # Documentation tools (optional)
    doxygen \
    graphviz \
    # Text editors for command-line work
    vim \
    nano \
    # Utilities for development
    tree \
    htop \
    lsof \
    strace \
    && rm -rf /var/lib/apt/lists/*

# Install uv CLI for Python tooling management
RUN curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 UV_INSTALL_DIR=/usr/local/bin sh

# Set up compiler alternatives to use clang-20/clang++-20/flang-20 by default
RUN update-alternatives --install /usr/bin/cc cc /usr/bin/clang-20 100 \
    && update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++-20 100 \
    && update-alternatives --install /usr/bin/clang clang /usr/bin/clang-20 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-20 100 \
    && update-alternatives --install /usr/bin/flang flang /usr/bin/flang-20 100 \
    && update-alternatives --install /usr/bin/lldb lldb /usr/bin/lldb-20 100 \
    && ln -sf /usr/bin/clang-20 /usr/bin/gcc \
    && ln -sf /usr/bin/clang++-20 /usr/bin/g++ \
    && ln -sf /usr/bin/flang-20 /usr/bin/gfortran

# Verify compiler versions
RUN clang --version \
    && clang++ --version \
    && flang --version \
    && cmake --version \
    && ninja --version

# User creation is handled by the common-utils devcontainer feature
# Set environment variables for the vscode user
ARG USERNAME=vscode
ENV PATH="/home/${USERNAME}/.local/bin:${PATH}"
ENV CC=clang-20
ENV CXX=clang++-20
ENV FC=flang-20

# Set the working directory
WORKDIR /workspaces

# User switching and git configuration will be handled by devcontainer features and postCreate script
# Default command
CMD ["/bin/bash"]
