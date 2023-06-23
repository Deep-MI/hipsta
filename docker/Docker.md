# Docker

Here we provide a [`Dockerfile`](Dockerfile) that can be used to create a Docker image and subsequently run the scripts within a Docker container. The main advantage of using virtualization technologies like Docker or Singularity is to create a controlled environment, with fixed versions of installable packages, which greatly helps with the stability and reproducibility of the analyses. Also, our Docker and Singularity images contain all additionally required external software packages (FreeSurfer and gmsh), which obviates the need for separate installations. Specifically, the Docker image that we provide here will be based on Ubuntu, contain the shape and thickness analysis scripts, a minimal subset of the FreeSurfer software, the gmsh software, and any additionally required packages from the Ubuntu distribution and the Python ecosystem.

## Build the Docker image

To build the Docker image, execute the following command after traversing into the *docker* directory of this repository:

```bash
docker build --rm -t hipsta -f Dockerfile .
```
The name of the image will be `hipsta`, and it will be built from the `Dockerfile` configuration file from the *docker* directory.

The `--rm` flag will remove intermediate containers after a successful build; `-t` specifies the name of the image, and `-f` indicates the configuration file from which to build.  You can also use `hipsta:<tag>` as the name of the image, with `<tag>` indicating an arbitrary version tag (default: latest).

Take a look at the contents of the [`Dockerfile`](Dockerfile) to see what is done during the build process: essentially, it is getting the Ubuntu 22.04 image, installing additional packages from the distribution, downloading the hipsta toolboxes, and setting the necessary environment variables. Unless the `Dockerfile` changes, the build process has to be done only once.

## Download the Docker image

As an alternative to building the Docker image yourself, you can also download our pre-built images from [Dockerhub](https://hub.docker.com/r/deepmi/hipsta/tags):

```bash
docker pull deepmi/hipsta:latest
```

## Run the Docker image

After building the Docker image, you can run it with the following command to see the help message of the main script:

```bash
docker run --rm --user XXXX:YYYY hipsta
```

* This corresponds to calling `python3 shapetools.py` from the command line for a non-dockerized version of the program.
* The `--rm` flag takes care of removing the container once the analysis finished.
* The `--user XXXX:YYYY` part should be changed to the appropriate user id (XXXX, a number) and group id (YYYY, also a number); both can be checked with the commands `id -u` and `id -g` on linux-like systems). All generated files will then belong to the specified user and group. Without the flag, the Docker container will be run as root with all corresponding privileges, which is strongly discouraged.
* You can run different versions of the image using `hipsta:<tag>` instead of `hipsta` and replacing `<tag>` with any particular version identifier.

Before running an analysis, you need to register at the FreeSurfer website (https://surfer.nmr.mgh.harvard.edu/registration.html) to acquire a valid license (for free). If you already have an existing FreeSurfer installation, you can use the license file of that installation as well.

An analysis can be performed by adding several options to the above command (and after adjusting the user-specific settings and file- and pathnames):

```bash
docker run \
    --rm \
    --user XXXX:YYYY \
    -v /path/to/my/input/directory:/path_to_input_directory_inside_docker \
    -v /path/to/my/output/directory:/path_to_output_directory_inside_docker \
    -v /path/to/my/freesurfer/installation/.license:/opt/freesurfer/.license    
    hipsta \
    --filename /path_to_input_directory_inside_docker/my_filename.mgz \
    --outputdir /path_to_output_directory_inside_docker \
    --hemi HEMISPHERE \
    --lut LOOKUP-TABLE \
    --logfiledir /path_to_output_directory_inside_docker
```

* The first two `-v` arguments mount your data directory and output directories into the Docker image (note that full, not relative, pathnames should be given). Inside the image, they are visible under the name following the colon (in this case `/path_to_input_directory_inside_docker` and `/path_to_output_directory_inside_docker`, but these can be different). From within the docker image / container, there will be read and write access to the directories that are mounted into the image (unless specified otherwise; for example, use `-v /path/to/my/input/directory:/path_to_input_directory_inside_docker:ro` for read-only access to your input directory).
* The third `-v` argument mounts your freesurfer license to the FreeSurfer directory within the Docker image (`/opt/freesurfer`). You can use a `.license` or a `license.txt` file, which can, but does not need to, be present in your local FreeSurfer directory on your system.
* The next part of the Docker command is the name of the Docker image, which is `hipsta`.
* After that, all other flags are identical to the ones that are used for the `shapetools.py` program (which are explained on the main page and the help message of the program). Note that, in addition to the `--filename`, `--hemi`, `--lut` arguments, which are mandatory, additional arguments such as `--outputdir` can be specified - in the same way as for non-dockerized version of the program. Also note that file- and pathnames need to correspond to the targets of the file / directory mappings within the Docker image, not to the local system. Finally, `HEMISPHERE` should either be `lh` or `rh`, and `LOOKUP-TABLE` should be `fs711` or `ashs`.
