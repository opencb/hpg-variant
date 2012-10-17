#!/bin/bash


if [ $# -lt 2 ]
then
	echo "local and remote branch required!!" 
	echo "./update-libs.sh next next-from-all"
	exit -1
fi

CURRENT_DIR=`pwd`

for i in libs/*;
do
	cd $i
	git checkout $1
	git pull origin $2
	cd $CURRENT_DIR
done
