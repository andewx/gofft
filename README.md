# fft [![GoDoc][godoc-badge]][godoc] [![Build Status][travis-ci-badge]][travis-ci] [![Report Card][report-card-badge]][report-card]

A better radix-2 fast Fourier transform in Go.

Package fft provides an efficient radix-2 fast discrete Fourier transformation algorithm in pure Go.

This code is much faster than existing FFT implementations and uses no additional memory.

The algorithm is non-recursive, works in-place overwriting the input array, and requires O(1) additional space.

## What

I took [an existing](https://github.com/ktye/fft) FFT implementation in Go, cleaned and improved the code API and performance, and replaced the permutation step with an algorithm that works with no temp array.

Performance was more than doubled over the original code, and is consistently the fastest Go FFT library (see benchmarks below) while remaining in pure Go.

Added convolution functions `Convolve(x, y)`, `FastConvolve(x, y)`, `MultiConvolve(x...)`, `FastMultiConvolve(X)`, which implement the discrete convolution and a new hierarchical convolution algorithm that has utility in a number of CS problems. This computes the convolution of many arrays in $O(n log^2 n)$ run time, and in the case of FastMultiConvolve O(1) additional space.

Also included new utility functions: `IsPow2`, `NextPow2`, `ZeroPad`, `ZeroPadToNextPow2`, `Float64ToComplex128Array`, `Complex128ToFloat64Array`, and `RoundFloat64Array`

## Why

Most existing FFT libraries in Go allocate temporary arrays with O(N) additional space. This is less-than-ideal when you have arrays of length of $2^{25}$ or more, where you quickly end up allocating gigabytes of data and dragging down the FFT calculation to a halt.

Additionally, the new convolution functions have significant utility for projects I've written or am planning.

One downside is that the FFT is not multithreaded (like go-dsp is), so for large vector size FFTs on a multi-core machine it will be slower than it could be. FFTs can be run in parallel, however, so in the case of many FFT calls it will be faster.

## How

```go
package main

import (
	"fmt"
	"github.com/argusdusty/fft"
)

func main() {
	// Do an FFT and IFFT and get the same result
	testArray := fft.Float64ToComplex128Array([]float64{1, 2, 3, 4, 5, 6, 7, 8})
	err := fft.FFT(testArray)
	if err != nil {
		panic(err)
	}
	err = fft.IFFT(testArray)
	if err != nil {
		panic(err)
	}
	result := fft.Complex128ToFloat64Array(testArray)
	fft.RoundFloat64Array(result)
	fmt.Println(result)

	// Do a discrete convolution of the testArray with itself
	testArray, err = fft.Convolve(testArray, testArray)
	if err != nil {
		panic(err)
	}
	result = fft.Complex128ToFloat64Array(testArray)
	fft.RoundFloat64Array(result)
	fmt.Println(result)
}
```

Outputs:

```text
[1 2 3 4 5 6 7 8]
[1 4 10 20 35 56 84 120 147 164 170 164 145 112 64]
```

### Benchmarks

```text
fft>go test -bench=FFT$ -benchmem -cpu=1 -benchtime=5s
goos: windows
goarch: amd64
pkg: github.com/argusdusty/fft
cpu: AMD Ryzen 9 5900X 12-Core Processor
```

| Algorithm      | Size               | Iterations | Time             | Throughput    | Memory         | Allocs       |
|:---------------|:-------------------|:-----------|:-----------------|:--------------|:---------------|:-------------|
| Naive          | Tiny (4)           | 43044590   | 149.1 ns/op      | 429.29 MB/s   | 64 B/op        | 1 allocs/op  |
| Naive          | Small (128)        | 31022      | 199949 ns/op     | 10.24 MB/s    | 2048 B/op      | 1 allocs/op  |
| Naive          | Medium (4096)      | 31         | 190195290 ns/op  | 0.34 MB/s     | 65536 B/op     | 1 allocs/op  |
| ktye           | Tiny (4)           | 334951483  | 16.13 ns/op      | 3968.29 MB/s  | 0 B/op         | 0 allocs/op  |
| ktye           | Small (128)        | 5627564    | 1064 ns/op       | 1923.91 MB/s  | 0 B/op         | 0 allocs/op  |
| ktye           | Medium (4096)      | 98692      | 62663 ns/op      | 1045.86 MB/s  | 0 B/op         | 0 allocs/op  |
| ktye           | Large (131072)     | 1436       | 4210853 ns/op    | 498.03 MB/s   | 0 B/op         | 0 allocs/op  |
| mjibson/go-dsp | Tiny (4)           | 1864420    | 2882 ns/op       | 22.21 MB/s    | 499 B/op       | 13 allocs/op |
| mjibson/go-dsp | Small (128)        | 701971     | 8934 ns/op       | 229.24 MB/s   | 5572 B/op      | 18 allocs/op |
| mjibson/go-dsp | Medium (4096)      | 33637      | 149844 ns/op     | 437.36 MB/s   | 164358 B/op    | 23 allocs/op |
| mjibson/go-dsp | Large (131072)     | 993        | 6568056 ns/op    | 319.30 MB/s   | 5243432 B/op   | 28 allocs/op |
| mjibson/go-dsp | Huge (4194304)     | 13         | 463206769 ns/op  | 144.88 MB/s   | 167772795 B/op | 33 allocs/op |
| mjibson/go-dsp | Massive (16777216) | 3          | 2478991133 ns/op | 108.28 MB/s   | 671089306 B/op | 35 allocs/op |
| gonum          | Tiny (4)           | 123918710  | 50.12 ns/op      | 1276.84 MB/s  | 0 B/op         | 0 allocs/op  |
| gonum          | Small (128)        | 2041134    | 2904 ns/op       | 705.33 MB/s   | 0 B/op         | 0 allocs/op  |
| gonum          | Medium (4096)      | 42642      | 138597 ns/op     | 472.85 MB/s   | 0 B/op         | 0 allocs/op  |
| gonum          | Large (131072)     | 870        | 6913295 ns/op    | 303.35 MB/s   | 0 B/op         | 0 allocs/op  |
| gonum          | Huge (4194304)     | 15         | 359428827 ns/op  | 186.71 MB/s   | 0 B/op         | 0 allocs/op  |
| gonum          | Massive (16777216) | 4          | 1407001550 ns/op | 190.79 MB/s   | 0 B/op         | 0 allocs/op  |
| scientificgo   | Tiny (4)           | 66733696   | 77.52 ns/op      | 825.58 MB/s   | 128 B/op       | 2 allocs/op  |
| scientificgo   | Small (128)        | 3675384    | 1653 ns/op       | 1238.64 MB/s  | 4096 B/op      | 2 allocs/op  |
| scientificgo   | Medium (4096)      | 88891      | 60442 ns/op      | 1084.28 MB/s  | 131072 B/op    | 2 allocs/op  |
| scientificgo   | Large (131072)     | 2526       | 2886025 ns/op    | 726.66 MB/s   | 4194304 B/op   | 2 allocs/op  |
| scientificgo   | Huge (4194304)     | 26         | 228963535 ns/op  | 293.10 MB/s   | 134217728 B/op | 2 allocs/op  |
| scientificgo   | Massive (16777216) | 4          | 1312125600 ns/op | 204.58 MB/s   | 536870914 B/op | 2 allocs/op  |
| fft          | Tiny (4)           | 1000000000 | 5.687 ns/op      | 11253.23 MB/s | 0 B/op         | 0 allocs/op  |
| fft          | Small (128)        | 6283482    | 1037 ns/op       | 1974.64 MB/s  | 0 B/op         | 0 allocs/op  |
| fft          | Medium (4096)      | 120750     | 51886 ns/op      | 1263.07 MB/s  | 0 B/op         | 0 allocs/op  |
| fft          | Large (131072)     | 2596       | 2387875 ns/op    | 878.25 MB/s   | 0 B/op         | 0 allocs/op  |
| fft          | Huge (4194304)     | 27         | 265444115 ns/op  | 252.82 MB/s   | 0 B/op         | 0 allocs/op  |
| fft          | Massive (16777216) | 5          | 1069123500 ns/op | 251.08 MB/s   | 0 B/op         | 0 allocs/op  |

[travis-ci-badge]:   https://api.travis-ci.org/argusdusty/fft.svg?branch=master
[travis-ci]:         https://api.travis-ci.org/argusdusty/fft
[godoc-badge]:       https://godoc.org/github.com/argusdusty/fft?status.svg
[godoc]:             https://godoc.org/github.com/argusdusty/fft
[report-card-badge]: https://goreportcard.com/badge/github.com/argusdusty/fft
[report-card]:       https://goreportcard.com/report/github.com/argusdusty/fft