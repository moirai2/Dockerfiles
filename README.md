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

## command lines
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
* To download GitHub
```
git pull https://github.com/moirai2/Dockerfiles.git
```
* To update GitHub
```
cd Dockerfiles
git pull
```
* To start docker-compose
```
docker-compose build
```
* To start web services
```
docker-compose up -d
```
* To stop web services
```
docker-compose down
```
* Remove orphans
```
docker-compose down --rmi all --volumes --remove-orphans
```