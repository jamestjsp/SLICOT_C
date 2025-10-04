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
    # Compilers (pinned to version 13 for consistency)
    gcc-13 \
    g++-13 \
    gfortran-13 \
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
    # Python for scripts and utilities
    python3 \
    python3-pip \
    python3-venv \
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

# Set up compiler alternatives to use gcc-13/g++-13/gfortran-13 by default
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-13 100 \
    && update-alternatives --install /usr/bin/gcov gcov /usr/bin/gcov-13 100

# Verify compiler versions
RUN gcc --version \
    && g++ --version \
    && gfortran --version \
    && cmake --version \
    && ninja --version

# User creation is handled by the common-utils devcontainer feature
# Set environment variables for the vscode user
ARG USERNAME=vscode
ENV PATH="/home/${USERNAME}/.local/bin:${PATH}"
ENV CC=gcc
ENV CXX=g++
ENV FC=gfortran

# Install Python packages (optional but useful for data analysis/visualization)
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    scipy \
    matplotlib \
    pandas

# Set the working directory
WORKDIR /workspaces

# User switching and git configuration will be handled by devcontainer features and postCreate script
# Default command
CMD ["/bin/bash"]
