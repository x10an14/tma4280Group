#Create CMAKE sub-folders
mkdir debug release -p

#Delete whatever was there before
rm -rf release/*

cd ../release/
CXX=icpc CC=icc FC=ifort cmake .. -DCMAKE_BUILD_TYPE=Release
make clean && make

cd ..
chmod 755 release/ex4
chmod 755 ex4_1n2ppnRel.sh
chmod 755 ex4_1n4ppnRel.sh
chmod 755 ex4_1n8ppnRel.sh
chmod 755 ex4_2n1ppnRel.sh
chmod 755 ex4_2n2ppnRel.sh
chmod 755 ex4_2n4ppnRel.sh

qsub ex4_1n2ppnRel.sh
qsub ex4_1n4ppnRel.sh
qsub ex4_1n8ppnRel.sh
qsub ex4_2n1ppnRel.sh
qsub ex4_2n2ppnRel.sh
qsub ex4_2n4ppnRel.sh
