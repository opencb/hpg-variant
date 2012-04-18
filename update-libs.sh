#!/bin/bash


for i in libs/*;
do
	cd $i;

	`git checkout $1`
	`git pull origin $1 `

	cd -;
done
