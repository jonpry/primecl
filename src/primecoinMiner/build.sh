#!/bin/sh
g++ clQueue.cpp oclMain.cpp bn2.cpp bn2_div.cpp -lgmpxx -lOpenCL -lcrypto -lssl -g
