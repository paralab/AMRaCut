set -e

mkdir -p build

cmake -G Ninja -S . -B build

ninja install -C ./build -v

echo -e "========== AMRaCut lib compile successfull ==========\n"



