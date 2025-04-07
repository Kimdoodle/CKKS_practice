git clone https://github.com/microsoft/SEAL.git
cd SEAL
cmake build -S . -B build
cmake --build build --config Debug
cmake --build build --config Release