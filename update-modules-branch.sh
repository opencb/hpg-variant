#!/bin/bash

BRANCH="next"
if [ $# -gt 0 ]; then
	BRANCH=$1
fi

for i in libs/*; do
	if [ -d "$i" ]; then
		cd $i
		git checkout $BRANCH
		cd -
	fi
done
