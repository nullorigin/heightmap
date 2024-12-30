# HighMap

A C++ library to generate two-dimensional terrain heightmaps for software rendering or video games. HighMap is the core solver of the [Hesiod](https://github.com/otto-link/Hesiod) GUI.

<p align="center"><img src="https://github.com/otto-link/HighMap/assets/121820229/4451f45b-c74e-4a06-9431-9302477a12c5" width="512"></p>

>[!WARNING] 
> Work in progress, use at your own risk!

## License

This project is licensed under the GNU General Public License v3.0.

## Getting started

### Getting the sources

Use `git` to retrieve the sources: 
``` bash
git clone git@github.com:otto-link/HighMap.git
cd HighMap
git submodule update --init --recursive
```

### Building

#### Linux

Build by making a build directory (i.e. `build/`), run `cmake` in that dir, and then use `make` to build the desired target.

Example:
``` bash
mkdir build && cd build
cmake ..
make
```

Simple examples are provided for most of the functions, for instance:
```
bin/./ex_noise_fbm
eog ex_noise_fbm_*.png
```

#### Windows

Use `PowerShell` to clone the repository (because Visual Studio won't clone the submodules) using `git` command lines (see above).

Install the missing OpenSource dependencies using `vcpkg`:
```
vcpkg install libpng glm ocl-icd-opencl-dev assimp opencv
```

You should then be able to build the sources using Visual Studio.

### Usage examples

See the `examples` folder.
