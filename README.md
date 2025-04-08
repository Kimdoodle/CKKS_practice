## 🪟 Installation Guide (Windows)

### ✅ Requirements
- Microsoft Visual Studio 2022

### ⚙️ Build Steps
1. Open ```x64 Native Tools Command Prompt for VS 2022```
  
2. Run the following commands:
```bash
git clone https://github.com/Kimdoodle/CKKS_practice.git SEAL_SGN
cd SEAL_SGN
git clone https://github.com/microsoft/SEAL.git
cd SEAL
cmake build -S . -B build
cmake --build build --config Release
```
3. Open ```SEAL_SGN``` folder in Visual Studio

4. Run ```SEAL_SGN.exe``` in x64 Release mode.
