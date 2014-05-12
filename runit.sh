#!/bin/bash
cd /home/justin/primecoin
while true
do
LD_PRELOAD="./libsleep.so" ./jhprimeminer -o http://ypool.net:10034 -u jonpry.foo -p x >> nohup.out
sleep 1
done
