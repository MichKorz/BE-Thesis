# Developer Cheatsheet: Environment Setup and LWE Compilation


## 0. Rebuilding the Container (if needed)

```bash
cd ~/Projects/BE-Thesis
podman build -t be-thesis-env .

```

---

## 1. Launching the Container Environment (Host)

To start the container with access to the RTX 5070 Ti GPU and mount the project directory with SELinux permissions (`:z`), run the following command in the Fedora system terminal:

```bash
podman run --rm -it \
  --device nvidia.com/gpu=all \
  -v ~/Projects/BE-Thesis:/workspace:z \
  be-thesis-env bash

```

---

## 2. Cloning / Updating Submodules (Host - Optional)

If the repository has been cloned on a new machine and the submodule folders are empty, initialize them on the host before entering the container:

```bash
git submodule update --init --recursive

```

---

## 3. Compiling Lightwave Explorer (Inside Container)

After launching the container, the terminal will automatically open in the `/workspace` directory. Execute the following commands:

```bash
# 1. Navigate to the LWE submodule directory
cd LightwaveExplorer

# 2. Create an isolated build directory and enter it
mkdir -p build && cd build

# 3. Clean previous configuration files (if any)
rm -rf *

# 4. Configure CMake for the Blackwell architecture (sm_120)[cite: 1]
cmake -DCLI=1 -DUSE_CUDA=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=120 ..

# 5. Compile the binaries using all available CPU cores[cite: 1]
cmake --build . -j$(nproc)

```

---

## 4. Preparing the Runtime Environment (Inside Container)

For the CLI executable to run properly, it requires the crystal database to be present in the same directory:

```bash
# Copy the required database to the build directory[cite: 1]
cp ../CrystalDatabase.txt .

```
