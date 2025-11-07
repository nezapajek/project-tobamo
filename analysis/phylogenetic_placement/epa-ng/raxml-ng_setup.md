
unzip raxml-ng_v1.2.2_linux_x86_64_MPI.zip -d raxml-ng_v1.2.2

cd raxml-ng_v1.2.2

# gcc ver 9.4.0
sudo apt install gcc 

# cmake ver 3.16.3
sudo apt install cmake

mkdir build && cd build

sudo apt install build-essential

cmake -DUSE_MPI=OFF ..

# once configuration succeeds
make -j 4

# test
./bin/raxml-ng --help

