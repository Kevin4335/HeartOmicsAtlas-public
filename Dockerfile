# Use R base image
FROM rocker/r-ver:4.3.1

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk-dev \
    libxt-dev \
    libx11-dev \
    libreadline-dev \
    zlib1g-dev \
    libicu-dev \
    libbz2-dev \
    liblzma-dev \
    pandoc \
    && apt-get clean

# Copy all project files
COPY . /app

# Install all R packages via script (fail build if anything goes wrong)
RUN Rscript install_packages.R

# Expose ports for each R server
EXPOSE 9025 9026 9027 9028

# Start all servers
CMD ["Rscript", "run_all.R"]
