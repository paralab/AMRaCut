set -e

mkdir -p build

cmake -G Ninja -S . -B build

ninja install -C ./build -v

echo -e "========== AMRaCut lib compile successfull ==========\n"

# cd test
# bash ./compile.sh

# mpirun -np 3 --oversubscribe ./build/test
# mpirun -np 1 --oversubscribe gdbserver :3000 ./build/test

# valgrind --tool=memcheck --leak-check=yes mpirun -np 1 ./build/test
# cmake --build ./build -v --target install

