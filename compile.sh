set -e

mkdir -p build

cmake -DAMRACUT_INTEGER_WIDTH=32 -S . -B build
cmake --build build --target install --verbose

# ====== if Ninja is installed ======

# cmake -G Ninja -S . -B build
# ninja install -C ./build -v

# ===================================

echo -e "========== AMRaCut lib compile successfull ==========\n"



