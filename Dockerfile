# Base image with Python, Julia, and R pre-installed
FROM registry.codeocean.com/codeocean/py-julia-r:python3.10.12-R4.2.3-julia1.7.0-ubuntu22.04

# Avoid prompts during the installation
ARG DEBIAN_FRONTEND=noninteractive
ARG MOSEKLM_LICENSE_FILE

# Install Python packages
RUN pip install -U --no-cache-dir \
    h5py==3.12.1 \
    matplotlib==3.10.0 \
    numpy==2.2.1 \
    pandas==2.2.3 \
    scipy==1.15.0

# Set MOSEK environment variables
ENV MOSEK_DIR=/opt/mosek
ENV MOSEKLM_LICENSE_FILE=/root/mosek/mosek.lic
ENV PATH="${MOSEK_DIR}/bin:${PATH}"

# Copy MOSEK license
RUN mkdir -p /root/mosek
COPY mosek.lic /root/mosek/mosek.lic

# Install Julia packages
RUN julia -e 'using Pkg; Pkg.add(PackageSpec(name="Combinatorics", version="1.0.2"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Dates"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="DynamicPolynomials", version="0.4.6"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="HDF5", version="0.17.2"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="JuMP", version="1.23.6"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="LinearAlgebra"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Mosek", version="10.2.0"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="MosekTools", version="0.15.1"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="NLopt", version="1.1.2"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="QuantumOptics", version="1.1.1"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Random"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Statistics"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/wangjie212/TSSOS"))' \
    && julia -e 'using Pkg; Pkg.API.precompile()'

# Copy files from GitHub repository
RUN git clone https://github.com/zpopovych/OQSID_for_Docker.git /root/OQSID
