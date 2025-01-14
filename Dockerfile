# Use a minimal base image with Ubuntu
FROM ubuntu:22.04

# Avoid prompts during installation
ARG DEBIAN_FRONTEND=noninteractive
ARG MOSEKLM_LICENSE_FILE

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    libffi-dev \
    libssl-dev \
    libreadline-dev \
    libbz2-dev \
    libsqlite3-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgdbm-dev \
    liblzma-dev \
    tk-dev \
    libxml2-dev \
    libxmlsec1-dev \
    xz-utils \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python 3.11.4
RUN wget https://www.python.org/ftp/python/3.11.4/Python-3.11.4.tgz \
    && tar -xvf Python-3.11.4.tgz \
    && cd Python-3.11.4 \
    && ./configure --enable-optimizations \
    && make -j$(nproc) \
    && make altinstall \
    && cd .. \
    && rm -rf Python-3.11.4 Python-3.11.4.tgz

# Install pip for Python
RUN python3.11 -m ensurepip --upgrade && python3.11 -m pip install --upgrade pip

# Install Julia 1.11.2
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.2-linux-x86_64.tar.gz \
    && tar -xvzf julia-1.11.2-linux-x86_64.tar.gz -C /usr/local --strip-components=1 \
    && rm julia-1.11.2-linux-x86_64.tar.gz

# Add Julia and Python to PATH
ENV PATH="/usr/local/bin:$PATH"

# Verify installations
RUN python3.11 --version && julia --version

# Set MOSEK environment variables
ENV MOSEK_DIR=/opt/mosek
ENV MOSEKLM_LICENSE_FILE=/root/mosek/mosek.lic
ENV PATH="${MOSEK_DIR}/bin:${PATH}"

# Copy MOSEK license
#RUN mkdir -p /root/mosek
#COPY mosek.lic /root/mosek/mosek.lic

# Install Python packages
RUN pip install -U --no-cache-dir \
    h5py==3.12.1 \
    matplotlib==3.10.0 \
    numpy==2.2.1 \
    pandas==2.2.3 \
    scipy==1.15.0

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
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="QuantumOptics"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Random"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(name="Statistics"))' \
    && julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/wangjie212/TSSOS"))' \
    && julia -e 'using Pkg; Pkg.API.precompile()'

# Copy files from GitHub repository
RUN git clone https://github.com/zpopovych/OQSID_for_Docker.git /root/OQSID

# Define a volume (optional, for documentation and persistence purposes)
VOLUME /root/OQSID/results

# Set the working directory
WORKDIR /root/OQSID/code

# Ensure the script is executable
RUN chmod +x run.sh

# Default command
CMD ["./run.sh"]