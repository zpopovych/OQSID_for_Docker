# Use the specified Miniconda base image
FROM registry.codeocean.com/codeocean/miniconda3:4.12.0-python3.9-ubuntu20.04

# Set non-interactive mode for apt-get
ARG DEBIAN_FRONTEND=noninteractive

# Define MOSEK environment variables
ENV MOSEK_DIR=/opt/mosek
ENV PATH="${MOSEK_DIR}/bin:${PATH}"

# Update and install required tools, Julia, Jupyter dependencies, and GUI libraries
RUN apt-get update || apt-get update --fix-missing \
    && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    build-essential \
    python3-pip \
    x11-apps \
    libx11-utils \
    xvfb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.3-linux-x86_64.tar.gz -O julia.tar.gz \
    && tar -xzf julia.tar.gz \
    && mv julia-1.9.3 /usr/local/julia \
    && ln -s /usr/local/julia/bin/julia /usr/local/bin/julia \
    && rm julia.tar.gz

# Download and install MOSEK
RUN wget https://download.mosek.com/stable/10.2.11/mosektoolslinux64x86.tar.bz2 -O mosek.tar.bz2 \
    && mkdir -p ${MOSEK_DIR} \
    && tar -xjf mosek.tar.bz2 -C ${MOSEK_DIR} --strip-components=1 \
    && rm mosek.tar.bz2

# Set up MOSEK license (ensure the license file is available in the build context)
COPY mosek.lic /root/.mosek/mosek.lic


# Set up MOSEK license (ensure the license file is available in the build context)
WORKDIR /code
COPY . /code

# Install Python Jupyter Notebook and required Python libraries
RUN pip install --no-cache-dir notebook h5py numpy matplotlib

# Install IJulia and other Julia packages for MOSEK and additional dependencies
RUN echo "using Pkg; \
    Pkg.add(\"IJulia\"); \
    Pkg.add(\"MosekTools\"); \
    Pkg.add(\"MathOptInterface\"); \
    Pkg.add(\"LinearAlgebra\"); \
    Pkg.add(\"QuantumOptics\"); \
    Pkg.add(\"DynamicPolynomials\"); \
    Pkg.add(\"Random\"); \
    Pkg.add(\"JuMP\"); \
    Pkg.add(\"NLopt\"); \
    Pkg.add(\"HDF5\"); \
    Pkg.add(\"Statistics\"); \
    Pkg.add(\"Dates\"); \
    Pkg.add(PackageSpec(url=\"https://github.com/wangjie212/TSSOS\"))" > setup.jl \
    && julia setup.jl \
    && rm setup.jl

# Expose Jupyter port
EXPOSE 8888

# Start Jupyter Notebook server with virtual display for GUI
CMD ["julia", "Hello.jl"]

