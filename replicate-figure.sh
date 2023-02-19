#!/bin/bash
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT
git clone https://github.com/hageldave/uamds.git
cd uamds/
git checkout replicability-stamp
cd java/core/
mvn install
cd ../plots/
mvn compile -P standalonejar assembly:single
mvn compile exec:java -D"exec.mainClass"="uamds.plots.Teaser" &
sleep 15s
cd ../../..
mv uamds/java/plots/teaser.svg teaser.svg

