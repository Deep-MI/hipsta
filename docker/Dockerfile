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
    python3-pip

# download FreeSurfer
RUN wget -qO- https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu22_x86_64-7.4.1.tar.gz | tar zxv --no-same-owner -C /opt \
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

# clone hipsta toolbox
RUN git clone https://github.com/Deep-MI/Hipsta.git /app/hipsta

# install additional python packages
RUN python3 -m pip install `cat /app/hipsta/requirements.txt`

# Add FreeSurfer environment variables (.license file needed, alternatively export FS_LICENSE=path/to/license)
ENV FREESURFER_HOME=/opt/freesurfer

# Add other environment variables
ENV OMP_NUM_THREADS=1

# Set the PATH
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/opt/gmsh-2.16.0-Linux/bin:${PATH}

# Set entrypoint (non-interactive mode)
ENTRYPOINT ["python3", "/app/hipsta/shapetools.py"]

## Run app when the container launches (interactive mode)
#CMD ["/bin/bash"]

