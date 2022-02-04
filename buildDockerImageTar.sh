name=$1
cd $name
docker build -t $name .
docker save $name | gzip -c > $name.tar.gz