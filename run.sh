#!/bin/sh
while true
do
	process=`ps -ef| grep FINSYS | grep -v grep | wc -l`
	if [ "$process" -eq 0 ];
	then
		echo "process down,now start"
		nohup ./FINSYS >/dev/null  2>log &
	else
		sleep 3;
	fi
done
