## 🪟 Installation Guide (Windows)

### ✅ Requirements
- Microsoft Visual Studio 2022

### ⚙️ Build Steps

1. In Visual Studio, create a CMake Project

 - Select ```Place solution and project in the same directory```

2. Remove ```CMakeLists.txt```, ```*.cpp, *.h``` in root directory
3. Open terminal and run the following commands:
```bash
git clone https://github.com/microsoft/SEAL.git
cd SEAL
cmake build -S . -B build
cmake --build build --config Release
git clone https://github.com/Kimdoodle/CKKS_practice.git .
```
4. Run ```SEAL_SGN.exe``` in x64 Release mode.