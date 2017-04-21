# blancmange
Compute [blancmange function](https://en.wikipedia.org/wiki/Blancmange_curve) of any positive rational number less than 1.

## Requirements
[GMP (GNU Multi-Precision Library)](https://gmplib.org/)

## Building
Install `libgmp` headers and libraries.

For example, on Ubuntu install packages:
```
libgmp10
libgmp-dev
```

Then compile using `make.sh`

## Usage

Program reads input number from `STDIN` and prints result to `STDOUT`.
Input number must be a positive decimal proper fraction (value must be between 0 and 1).

Example:
```
$ echo '1/5' | blancmange
8/15
```
