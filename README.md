# Primal-dual Interior-point Solver for Linear Programming

This is a course project for CS520: Computational Methods for Optimization collaborated with [Xin Ye](https://github.com/xinye83). The complete algorithms and experimental result are in the [report](./report.pdf).

## Features

- Julia implementation

- Presolve stage

- Mehrotra's predictor-corrector algorithm

## Requirements

- [Julia](https://julialang.org/) (v1.1)

- [MatrixDepot](https://github.com/JuliaMatrices/MatrixDepot.jl) (v0.8): test matrix collection for Julia

## Usage

``` shell
git clone https://github.com/dlguo/primal-dual-interior-point.git
cd primal-dual-interior-point
julia test.jl
```

In the interactive interface, choose 9 preset test problem in `LPnetlib` or use your own dataset.
