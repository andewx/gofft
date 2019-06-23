# gofft [![GoDoc][godoc-badge]][godoc] [![Build Status][travis-ci-badge]][travis-ci] [![Report Card][report-card-badge]][report-card]
A better radix-2 fast Fourier transform in Go.

Package gofft provides an efficient radix-2 fast discrete Fourier transformation algorithm in pure Go.

This code is much faster than existing FFT implementations and uses no additional memory.

The algorithm is non-recursive, works in-place overwriting the input array, and requires O(1) additional space.

## What
I took [an existing](https://github.com/ktye/fft) FFT implementation in Go, cleaned and improved the code API and performance, and replaced the permutation step with an algorithm that works with no temp array.

Performance was nearly doubled over the original code, and is many times faster than go-dsp for small inputs (see benchmarks in fft_test).

Added convolution functions `Convolve(x, y)`, `FastConvolve(x, y)`, `MultiConvolve(x...)`, `FastMultiConvolve(X)`, which implement the discrete convolution and a new hierarchical convolution algorithm that has utility in a number of CS problems. This computes the convolution of many arrays in O(n\*ln(n)<sup>2</sup>) run time, and in the case of FastMultiConvolve O(1) additional space.

Also included a new utility functions: `ZeroPad(x, N)`, `Float64ToComplex128Array`, and `Complex128ToFloat64Array` which are self-documenting.

## Why
Most existing FFT libraries in Go allocate temporary arrays with O(N) additional space. This is less-than-ideal when you have arrays of length of 2<sup>25</sup> or more, where you quickly end up allocating gigabytes of data and dragging down the FFT calculation to a halt.

This code is much faster than major existing fft libraries while reducing memory usage and remaining in pure Go.

Additionally, the new convolution functions have significant utility for projects I've written or am planning.

One downside is that the FFT is not multithreaded (like go-dsp is), so for large vector size FFTs on a multi-core machine it will be slower than it could be. FFTs can be run in parallel, however, so in the case of many FFT calls it will be faster.

## How
```go
package main

import (
	"fmt"
	"github.com/argusdusty/gofft"
)

func main() {
	// Do an FFT and IFFT and get the same result
	testArray := gofft.Float64ToComplex128Array([]float64{1, 2, 3, 4, 5, 6, 7, 8})
	gofft.Prepare(len(testArray))
	err := gofft.FFT(testArray)
	if err != nil {
		panic(err)
	}
	err = gofft.IFFT(testArray)
	if err != nil {
		panic(err)
	}
	result := gofft.Complex128ToFloat64Array(testArray)
	gofft.RoundFloat64Array(result)
	fmt.Println(result)

	// Do a discrete convolution of the testArray with itself
	testArray, err = gofft.Convolve(testArray, testArray)
	if err != nil {
		panic(err)
	}
	result = gofft.Complex128ToFloat64Array(testArray)
	gofft.RoundFloat64Array(result)
	fmt.Println(result)
}
```

Outputs:
```
[1 2 3 4 5 6 7 8]
[1 4 10 20 35 56 84 120 147 164 170 164 145 112 64]
```

### Benchmarks
```
github.com\argusdusty\gofft>go test -bench=. -benchmem -cpu=1,4
goos: windows
goarch: amd64
pkg: github.com/argusdusty/gofft
BenchmarkSlowFFT/Tiny_(4)                5000000               259 ns/op         246.56 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT/Tiny_(4)-4             10000000               249 ns/op         256.19 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT/Small_(128)                3000            355391 ns/op           5.76 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT/Small_(128)-4              5000            355787 ns/op           5.76 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT/Medium_(4096)                 3         349252400 ns/op           0.19 MB/s       65536 B/op          1 allocs/op
BenchmarkSlowFFT/Medium_(4096)-4               3         347247300 ns/op           0.19 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)               10000000               131 ns/op         486.48 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)-4             20000000               116 ns/op         551.20 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)              300000              4297 ns/op         476.56 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)-4            300000              4214 ns/op         485.94 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)             10000            180892 ns/op         362.29 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)-4           10000            171514 ns/op         382.10 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)              100          10477746 ns/op         200.15 MB/s     2097152 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)-4            100          10248000 ns/op         204.64 MB/s     2097154 B/op          1 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)                300000              3497 ns/op          18.30 MB/s         568 B/op         13 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)-4              500000              3408 ns/op          18.78 MB/s         504 B/op         14 allocs/op
BenchmarkGoDSPFFT/Small_(128)             200000             10960 ns/op         186.85 MB/s        5600 B/op         18 allocs/op
BenchmarkGoDSPFFT/Small_(128)-4           100000             19320 ns/op         106.00 MB/s        5785 B/op         34 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)            10000            170892 ns/op         383.49 MB/s      164380 B/op         23 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)-4          10000            195829 ns/op         334.66 MB/s      164831 B/op         54 allocs/op
BenchmarkGoDSPFFT/Large_(131072)             200           7961076 ns/op         263.43 MB/s     5243448 B/op         28 allocs/op
BenchmarkGoDSPFFT/Large_(131072)-4           300           5679089 ns/op         369.28 MB/s     5244217 B/op         74 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)               2         719410200 ns/op          93.28 MB/s   167772816 B/op         33 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)-4             3         371179166 ns/op         180.80 MB/s   167774186 B/op         98 allocs/op
BenchmarkFFT/Tiny_(4)                   30000000              42.1 ns/op        1521.16 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Tiny_(4)-4                 30000000              42.1 ns/op        1520.90 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)                 1000000              1304 ns/op        1570.10 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)-4               1000000              1273 ns/op        1608.04 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)                 20000             70867 ns/op         924.77 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)-4               20000             70953 ns/op         923.65 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)                  300           4475215 ns/op         468.61 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)-4                300           4474826 ns/op         468.66 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)                    3         453284966 ns/op         148.05 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)-4                  3         456955133 ns/op         146.86 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Tiny_(4)           30000000              38.2 ns/op        1676.16 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Tiny_(4)-4        100000000              13.4 ns/op        4767.56 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Small_(128)         1000000              1307 ns/op        1565.95 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Small_(128)-4       3000000               478 ns/op        4281.00 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Medium_(4096)         20000             73033 ns/op         897.34 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Medium_(4096)-4       50000             24109 ns/op        2718.32 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Large_(131072)          300           4449412 ns/op         471.33 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Large_(131072)-4       1000           1914525 ns/op        1095.39 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Huge_(4194304)            3         447525700 ns/op         149.96 MB/s          37 B/op          1 allocs/op
BenchmarkFFTParallel/Huge_(4194304)-4          5         270185780 ns/op         248.38 MB/s          41 B/op          1 allocs/op
```

[travis-ci-badge]:   https://api.travis-ci.org/argusdusty/gofft.svg?branch=master
[travis-ci]:         https://api.travis-ci.org/argusdusty/gofft
[godoc-badge]:       https://godoc.org/github.com/argusdusty/gofft?status.svg
[godoc]:             https://godoc.org/github.com/argusdusty/gofft
[report-card-badge]: https://goreportcard.com/badge/github.com/argusdusty/gofft
[report-card]:       https://goreportcard.com/report/github.com/argusdusty/gofft