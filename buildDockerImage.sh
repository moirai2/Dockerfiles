name=$1
cd $name
docker build -t $name . > stdout.txt 2> stderr.txt
