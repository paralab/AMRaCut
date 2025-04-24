set -e

mkdir -p build

cmake -G Ninja -S . -B build

ninja -C ./build

echo "======== octree-test compile successfull =========="

mpirun -np 4 --oversubscribe ./build/test-octree