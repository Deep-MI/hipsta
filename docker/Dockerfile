FROM ubuntu:22.04

# update OS
RUN apt-get update

# get additional packages
RUN apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    tar \
    zip \
    tcsh \
    git \
    libglu1-mesa-dev \
    libxrender-dev \
    libxcursor-dev \
    libxft-dev \
    libxinerama-dev \
    libgomp1 \
    libsm-dev \
    python3 \
    python3-pip \
    libnss3 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libcups2 \
    libdrm2 \
    libxkbcommon0 \
    libxcomposite1 \
    libxdamage1 \
    libxfixes3 \
    libxrandr2 \
    libgbm1 \
    libasound2 \
    libpango-1.0-0 \
    libcairo2

# download FreeSurfer
#RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu22_x86_64-7.4.1.tar.gz | tar zxv --no-same-owner -C /opt \
#    --exclude='freesurfer/trctrain' \
#    --exclude='freesurfer/subjects/fsaverage_sym' \
#    --exclude='freesurfer/subjects/fsaverage3' \
#    --exclude='freesurfer/subjects/fsaverage4' \
#    --exclude='freesurfer/subjects/fsaverage5' \
#    --exclude='freesurfer/subjects/fsaverage6' \
#    --exclude='freesurfer/subjects/cvs_avg35' \
#    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
#    --exclude='freesurfer/subjects/V1_average' \
#    --exclude='freesurfer/average/mult-comp-cor'
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz | tar zxv --no-same-owner -C /opt \
    --exclude='freesurfer/trctrain' \
    --exclude='freesurfer/subjects/fsaverage_sym' \
    --exclude='freesurfer/subjects/fsaverage3' \
    --exclude='freesurfer/subjects/fsaverage4' \
    --exclude='freesurfer/subjects/fsaverage5' \
    --exclude='freesurfer/subjects/fsaverage6' \
    --exclude='freesurfer/subjects/cvs_avg35' \
    --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
    --exclude='freesurfer/subjects/V1_average' \
    --exclude='freesurfer/average/mult-comp-cor'

# download gmsh
RUN wget -qO- https://gmsh.info/bin/Linux/gmsh-2.16.0-Linux64.tgz | tar zxv --no-same-owner -C /opt

#
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install git+https://github.com/Deep-MI/Hipsta.git

# Install Chrome for Kaleido/plotly QC image rendering. Installed to a
# fixed, world-readable location (rather than the default $HOME) and
# pointed to via BROWSER_PATH, since the image is intended to be run
# with an arbitrary --user; HOME is also set to a writable default for
# the same reason (Chrome needs a writable profile/user-data directory).
RUN mkdir -p /opt/chrome-for-testing \
    && plotly_get_chrome -y --path /opt/chrome-for-testing \
    && chmod -R a+rX /opt/chrome-for-testing
ENV BROWSER_PATH=/opt/chrome-for-testing/chrome-linux64/chrome
ENV HOME=/tmp

# Add FreeSurfer environment variables (.license file needed, alternatively export FS_LICENSE=path/to/license)
ENV FREESURFER_HOME=/opt/freesurfer

# Add other environment variables
ENV OMP_NUM_THREADS=1

# Set the PATH
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/opt/gmsh-2.16.0-Linux/bin:${PATH}

# Set entrypoint (non-interactive mode)
ENTRYPOINT ["run_hipsta"]

## Run app when the container launches (interactive mode)
#CMD ["/bin/bash"]

