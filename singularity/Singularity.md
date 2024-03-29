# Singularity

We host releases of the hipsta package and all required external packages as Docker images on [Dockerhub](https://hub.docker.com/r/deepmi/hipsta/tags). For use on HPC systems or in other cases where Docker is not preferred you can easily create a Singularity image from the Docker images. 

## Singularity Image Creation

For creating a singularity image from the Dockerhub just run: 

```bash
cd /home/user/my_singlarity_images
singularity build hipsta-latest.sif docker://deepmi/hipsta:latest
```

Singularity Images are saved as `.sif` files. Here the _/homer/user/my_singlarity_images_ is the path where you want your file saved. You can change _deepmi/hipsta:latest_ with any tag provided in our [Dockerhub](https://hub.docker.com/r/deepmi/hipsta/tags).

If you want to use a locally available image that you created yourself, instead run:

```bash
cd /home/user/my_singlarity_images
singularity build hipsta-myimage.sif docker-daemon://hipsta:myimage
```

For how to create your own Docker images see our [Docker guide](../docker/Docker.md)

## Singularity Image Usage

After building the Singularity image, you need to register at the FreeSurfer website (https://surfer.nmr.mgh.harvard.edu/registration.html) to acquire a valid license (for free). If you already have an existing FreeSurfer installation, you can use the license file of that installation as well.

To run the shape and thickness analysis on a given subject using the Singularity image, execute the following command (after adjusting to your particular environment):

```bash
singularity exec \
    -B /path/to/my/input/directory:/path_to_input_directory_inside_container \
    -B /path/to/my/output/directory:/path_to_output_directory_inside_container \
    -B /path/to/my/freesurfer/installation/.license:/opt/freesurfer/.license
    /home/user/my_singularity_images/hipsta.sif \
    run_hipsta \
    --filename /path_to_input_directory_inside_container/my_filename.mgz \
    --outputdir /path_to_output_directory_inside_container \
    --hemi HEMISPHERE \
    --lut LOOKUP-TABLE
```

* The first two `-B` arguments mount your data directory and output directories into the singularity image (note that full, not relative, pathnames should be given). Inside the image, they are visible under the name following the colon (in this case `/path_to_input_directory_inside_container` and `/path_to_output_directory_inside_container`, but these can be different). From within the singularity image / container, there will be read and write access to the directories that are mounted into the image (unless specified otherwise).
* The third `-B` argument mounts your freesurfer license to the FreeSurfer directory within the singularity image (`/opt/freesurfer`). You can use a `.license` or a `license.txt` file, which can, but does not need to, be present in your local FreeSurfer directory on your system. The version of the license should correspond to the version inside the [`Dockerfile`](../docker/Dockerfile) that was used to create the original Docker image.
* The next part of the command is the name of the Singularity image, which is `hipsta.sif`. In this example it is located in `/home/user/my_singularity_images`, but the specific path will likely be different on your local system.
* For the Singularity image, we also have to explicitly specify the command that we want run, i.e. `run_hipsta`.
* After that, all other flags are identical to the ones that are used for the `run_hipsta` program (which are explained on the main page and the help message of the program). Note that, in addition to the `--filename`, `--hemi`, `--lut` arguments, which are mandatory, additional arguments such as `--outputdir` can be specified. Also note that file- and pathnames need to correspond to the targets of the file / directory mappings within the singularity image, not to the local system. Finally, `HEMISPHERE` should either be `lh` or `rh`, and `LOOKUP-TABLE` should be `freesurfer`, `ashs-penn_abc_3t_t2`, `ashs-umcutrecht_7t` or a custom segmentation.
