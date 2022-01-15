# Dockerfiles
* Collection of Dockerfiles I used for moirai2 project.
* Command to build docker image:
```
cd DIR
docker build -t NAME .
```
* Dockerhub: https://hub.docker.com
  * username: moirai2
  * moirai2/skewc: https://hub.docker.com/repository/docker/moirai2/skewc

## binaries
* Creates a Docker image
```
cd Dockerfiles
buildDockerImage.sh IMGNAME
```
* Creates a Docker image tar
```
cd Dockerfiles
buildDockerImageTar.sh IMGNAME
```
* Creates a singularity image file (SIF)
```
cd Dockerfiles
buildSingularityImageFile.sh IMGNAME
```