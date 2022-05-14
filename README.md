# GenCut

## Introduction

This is a rewrite of the [ScerSeg](https://github.com/mgeorgoulopoulos/ScerSeg) algorithm, which allows to segment a genome based on the 3-dimensional positions of genes, in such a way that the resulting segments (or clusters) form groups of statistical significance against a "signal", which is a property assigned to each gene.

The rewrite:

* Is in plain C++. No Qt, no external libraries required, except rapidcsv, which is automatically downloaded and built by the cmake script.
* Is a lot user-friendlier. Takes csv files as input instead of a SQLite database.
* Abstracts the signal source, allowing you to use lists (one scalar value assigned to each gene) and matrices (numeric value assigned to each pair of genes).

## Building

Prerequisites

* You need a build environment for C++, for which cmake can generate project files. This typically means Visual Studio (windows) and make/gcc (linux).
* Cmake must be installed
* Git must be installed (to fetch rapidcsv).

### GCC build

```
 mkdir bin
 cd bin
 cmake ../src
 make
 ./gencut --settings ../SampleSettings.csv
```

### Windows build

* Run cmake GUI and point it to src directory. Configure and generate.
* Open generated VS solution and build.

