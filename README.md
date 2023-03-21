# GenAlg
Generalised genetic algorithm library. Written in C with sacrifices to code-cleaniness for targetting memory and CPU efficiency.

# OS Support
Can be built with CMake. Tested and working on Windows, Linux, and FreeBSD.

# How To Build:
```
mkdir build && cd build
cmake (OPTIONAL ARCH FLAG: e.g. "-A x64") ..
cmake --build . (ADD: --target [genalg test usage_example] to build a subset) (ADD: "--config Release" for Windows build)
```
