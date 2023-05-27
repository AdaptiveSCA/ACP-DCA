# ACP-DCA against SE-SPECK64-K128

This is a script for ACP-DCA against a self-equivalence white-box SPECK implementation. The SE-SPECK is constructed with block size 64 and key size 128. The encodings are affine self-equivalences.

# Experiment Environment
Windows 10

Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz   3.70 GHz

16.0 GB RAM

gcc version 8.2.0 (MinGW.org GCC-8.2.0-3)

cmake version 3.18.3

# Build

```
$ gcc default_white_box_speck_64_K128.c
```

## Run

```
$ .\a.exe
```