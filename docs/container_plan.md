# Container image status and plan

The Mindboggle container (`nipy/mindboggle`) still runs locally with `mindboggle -h`, so there's still a working container that can be used as a baseline for running the full pipeline today, even though the image is old and large.

Rebuilding the container from `install/neurodocker.sh` is currently broken. The script is based on `neurodebian:stretch`, and the generated Dockerfile fails during `apt-get update` because the Debian Stretch repositories it depends on are no longer available on the normal mirrors. This means the existing published image is usable, but the source build path can't be reproduced as it is right now. 

The current image is also much larger than it needs to be for the long term. It bundles FreeSurfer, ANTs, Mindboggle, templates, and a large scientific Python stack into one image, which pushes the total size above 6 GB. Cleaning conda and apt caches, and removing cloned Git metadata after installation could help with the size. Because niwrap already integrates ANTs and FreeSurfer, though, when we integrate with niwrap, the image size will be reduced significantly. 

The first priority for modernization is to move the container build off Debian Stretch and onto a supported base image. Once the build works again, the next step should be to apply the safe cleanup reductions above and confirm that the image still runs correctly. After that, we can move towards a slimmer niwrap-based setup.