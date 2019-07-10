# gofft [![GoDoc][godoc-badge]][godoc] [![Build Status][travis-ci-badge]][travis-ci] [![Report Card][report-card-badge]][report-card]
A better radix-2 fast Fourier transform in Go.

Package gofft provides an efficient radix-2 fast discrete Fourier transformation algorithm in pure Go.

This code is much faster than existing FFT implementations and uses no additional memory.

The algorithm is non-recursive, works in-place overwriting the input array, and requires O(1) additional space.

## What
I took [an existing](https://github.com/ktye/fft) FFT implementation in Go, cleaned and improved the code API and performance, and replaced the permutation step with an algorithm that works with no temp array.

Performance was more than doubled over the original code, and is consistently the fastest Go FFT library (see benchmarks below)

Added convolution functions `Convolve(x, y)`, `FastConvolve(x, y)`, `MultiConvolve(x...)`, `FastMultiConvolve(X)`, which implement the discrete convolution and a new hierarchical convolution algorithm that has utility in a number of CS problems. This computes the convolution of many arrays in O(n\*ln(n)<sup>2</sup>) run time, and in the case of FastMultiConvolve O(1) additional space.

Also included new utility functions: `IsPow2`, `NextPow2`, `ZeroPad`, `ZeroPadToNextPow2`, `Float64ToComplex128Array`, `Complex128ToFloat64Array`, and `RoundFloat64Array`

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
BenchmarkSlowFFT/Tiny_(4)                       10000000               245 ns/op         261.18 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT/Tiny_(4)-4                     10000000               243 ns/op         263.11 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT/Small_(128)                        5000            320935 ns/op           6.38 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT/Small_(128)-4                      5000            319540 ns/op           6.41 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT/Medium_(4096)                         5         313555940 ns/op           0.21 MB/s       65536 B/op          1 allocs/op
BenchmarkSlowFFT/Medium_(4096)-4                       5         314757120 ns/op           0.21 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)                       20000000               119 ns/op         534.54 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)-4                     20000000               105 ns/op         608.55 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)                      500000              3743 ns/op         547.02 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)-4                    500000              3881 ns/op         527.62 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)                     10000            175428 ns/op         373.58 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)-4                   10000            165554 ns/op         395.86 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)                      200           8886224 ns/op         236.00 MB/s     2097152 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)-4                    200           9300254 ns/op         225.49 MB/s     2097154 B/op          1 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)                        500000              3219 ns/op          19.88 MB/s         545 B/op         13 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)-4                      500000              3791 ns/op          16.88 MB/s         504 B/op         14 allocs/op
BenchmarkGoDSPFFT/Small_(128)                     200000             12431 ns/op         164.74 MB/s        5598 B/op         18 allocs/op
BenchmarkGoDSPFFT/Small_(128)-4                   100000             22410 ns/op          91.39 MB/s        5785 B/op         34 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)                    10000            182113 ns/op         359.86 MB/s      164383 B/op         23 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)-4                  10000            208641 ns/op         314.11 MB/s      164830 B/op         54 allocs/op
BenchmarkGoDSPFFT/Large_(131072)                     200           8108197 ns/op         258.65 MB/s     5243448 B/op         28 allocs/op
BenchmarkGoDSPFFT/Large_(131072)-4                   300           6190119 ns/op         338.79 MB/s     5244211 B/op         74 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)                       2         893609500 ns/op          75.10 MB/s    167772840 B/op        33 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)-4                     3         399930300 ns/op         167.80 MB/s    167774058 B/op        97 allocs/op
BenchmarkScientificFFT/Tiny_(4)                 10000000               120 ns/op         532.54 MB/s         128 B/op          2 allocs/op
BenchmarkScientificFFT/Tiny_(4)-4               20000000               115 ns/op         555.60 MB/s         128 B/op          2 allocs/op
BenchmarkScientificFFT/Small_(128)               1000000              2278 ns/op         898.68 MB/s        4096 B/op          2 allocs/op
BenchmarkScientificFFT/Small_(128)-4             1000000              2172 ns/op         942.84 MB/s        4096 B/op          2 allocs/op
BenchmarkScientificFFT/Medium_(4096)               20000             88960 ns/op         736.69 MB/s      131072 B/op          2 allocs/op
BenchmarkScientificFFT/Medium_(4096)-4             20000             84723 ns/op         773.53 MB/s      131072 B/op          2 allocs/op
BenchmarkScientificFFT/Large_(131072)                300           4451421 ns/op         471.12 MB/s     4194304 B/op          2 allocs/op
BenchmarkScientificFFT/Large_(131072)-4              300           4597706 ns/op         456.13 MB/s     4194304 B/op          2 allocs/op
BenchmarkScientificFFT/Huge_(4194304)                  5         232384120 ns/op         288.78 MB/s    134217728 B/op         2 allocs/op
BenchmarkScientificFFT/Huge_(4194304)-4                5         223402620 ns/op         300.39 MB/s    134217728 B/op         2 allocs/op
BenchmarkFFT/Tiny_(4)                           30000000              44.9 ns/op        1424.97 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Tiny_(4)-4                         30000000              45.6 ns/op        1404.16 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)                         1000000              1493 ns/op        1371.70 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)-4                       1000000              1495 ns/op        1368.99 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)                         20000             83526 ns/op         784.61 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)-4                       20000             84523 ns/op         775.36 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)                          300           3972706 ns/op         527.89 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)-4                        500           4328479 ns/op         484.50 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)                            5         230782600 ns/op         290.79 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)-4                          5         230987220 ns/op         290.53 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Tiny_(4)                   30000000              44.8 ns/op        1429.20 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Tiny_(4)-4                 50000000              36.9 ns/op        1736.23 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Small_(128)                 1000000              1497 ns/op        1367.16 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Small_(128)-4               3000000               522 ns/op        3921.35 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Medium_(4096)                 20000             83676 ns/op         783.21 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Medium_(4096)-4               50000             26788 ns/op        2446.44 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Large_(131072)                  500           3967385 ns/op         528.60 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Large_(131072)-4               1000           1439149 ns/op        1457.22 MB/s           0 B/op          0 allocs/op
BenchmarkFFTParallel/Huge_(4194304)                    5         216819940 ns/op         309.51 MB/s          35 B/op          1 allocs/op
BenchmarkFFTParallel/Huge_(4194304)-4                 10         159174100 ns/op         421.61 MB/s          20 B/op          0 allocs/op
```

[travis-ci-badge]:   https://api.travis-ci.org/argusdusty/gofft.svg?branch=master
[travis-ci]:         https://api.travis-ci.org/argusdusty/gofft
[godoc-badge]:       https://godoc.org/github.com/argusdusty/gofft?status.svg
[godoc]:             https://godoc.org/github.com/argusdusty/gofft
[report-card-badge]: https://goreportcard.com/badge/github.com/argusdusty/gofft
[report-card]:       https://goreportcard.com/report/github.com/argusdusty/gofft