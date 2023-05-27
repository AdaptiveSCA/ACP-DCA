# ACP-DCA against IF-SPECK64-K128

This is a script for ACP-DCA against an implicit white-box SPECK implementation. The IF-SPECK is constructed with block size 64 and key size 128. The encodings are affine-quadratic self-equivalences. The degree of the implicit function is 3.

# Experiment Environment
Windows 10

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

cmake version 3.18.3

# Build

```
$ mkdir build
$ cd build
$ cmake -DBUILD_TARGET=win32 -G "MinGW Makefiles" ..
$ mingw32-make
```

## Run

```
$ .\ACPDCA.exe
```