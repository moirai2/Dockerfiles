name=$1
cd $name
docker build -t $name .
