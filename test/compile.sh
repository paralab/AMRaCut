set -e

mkdir -p build

cmake -G Ninja -S . -B build

ninja -C ./build

echo "======== test driver compile successfull =========="

mpirun -np 3 --oversubscribe ./build/test