# gofft [![GoDoc][godoc-badge]][godoc] [![Build Status][travis-ci-badge]][travis-ci] [![Report Card][report-card-badge]][report-card]
A better radix-2 fast Fourier transform in Go.

Package gofft provides an efficient radix-2 fast discrete Fourier transformation algorithm in pure Go.

This code is much faster than existing FFT implementations and uses no additional memory.

The algorithm is non-recursive, works in-place overwriting the input array, and requires O(1) additional space.

## What
I took [an existing](https://github.com/ktye/fft) FFT implementation in Go, cleaned and improved the code API and performance, and replaced the permutation step with an algorithm that works with no temp array.

Performance was more than doubled over the original code, and is consistently the fastest Go FFT library (see benchmarks below) while remaining in pure Go.

Added convolution functions `Convolve(x, y)`, `FastConvolve(x, y)`, `MultiConvolve(x...)`, `FastMultiConvolve(X)`, which implement the discrete convolution and a new hierarchical convolution algorithm that has utility in a number of CS problems. This computes the convolution of many arrays in O(n\*ln(n)<sup>2</sup>) run time, and in the case of FastMultiConvolve O(1) additional space.

Also included new utility functions: `IsPow2`, `NextPow2`, `ZeroPad`, `ZeroPadToNextPow2`, `Float64ToComplex128Array`, `Complex128ToFloat64Array`, and `RoundFloat64Array`

## Why
Most existing FFT libraries in Go allocate temporary arrays with O(N) additional space. This is less-than-ideal when you have arrays of length of 2<sup>25</sup> or more, where you quickly end up allocating gigabytes of data and dragging down the FFT calculation to a halt.

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
gofft>go test -bench=FFT$ -benchmem -cpu=1 -benchtime=5s
goos: windows
goarch: amd64
pkg: github.com/argusdusty/gofft
BenchmarkSlowFFT/Tiny_(4)                     25106818               233 ns/op         274.49 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT/Small_(128)                     18020            331225 ns/op           6.18 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT/Medium_(4096)                      19         325186379 ns/op           0.20 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)                     99976838              64.3 ns/op         995.83 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)                   2542816              2359 ns/op         868.34 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)                   55905            110228 ns/op         594.55 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)                    651           9742529 ns/op         215.26 MB/s     2097152 B/op          1 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)                     1787180              3383 ns/op          18.92 MB/s         522 B/op         13 allocs/op
BenchmarkGoDSPFFT/Small_(128)                   537966             11952 ns/op         171.36 MB/s        5589 B/op         18 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)                  34041            167400 ns/op         391.49 MB/s      164373 B/op         23 allocs/op
BenchmarkGoDSPFFT/Large_(131072)                   885           6844395 ns/op         306.40 MB/s     5243448 B/op         28 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)                     6         934443067 ns/op          71.82 MB/s   167772810 B/op         33 allocs/op
BenchmarkGoDSPFFT/Massive_(16777216)                 2        3758176900 ns/op          71.43 MB/s   671089328 B/op         35 allocs/op
BenchmarkGonumFFT/Tiny_(4)                    80345838              75.2 ns/op         851.17 MB/s           0 B/op          0 allocs/op
BenchmarkGonumFFT/Small_(128)                  1256538              4732 ns/op         432.75 MB/s           0 B/op          0 allocs/op
BenchmarkGonumFFT/Medium_(4096)                  25922            233528 ns/op         280.63 MB/s           0 B/op          0 allocs/op
BenchmarkGonumFFT/Large_(131072)                   526          11581403 ns/op         181.08 MB/s           0 B/op          0 allocs/op
BenchmarkGonumFFT/Huge_(4194304)                    12         485502825 ns/op         138.23 MB/s           0 B/op          0 allocs/op
BenchmarkGonumFFT/Massive_(16777216)                 3        2137556667 ns/op         125.58 MB/s           0 B/op          0 allocs/op
BenchmarkScientificFFT/Tiny_(4)               52617824               118 ns/op         542.92 MB/s         128 B/op          2 allocs/op
BenchmarkScientificFFT/Small_(128)             2982244              1896 ns/op        1080.37 MB/s        4096 B/op          2 allocs/op
BenchmarkScientificFFT/Medium_(4096)             90750             69664 ns/op         940.75 MB/s      131072 B/op          2 allocs/op
BenchmarkScientificFFT/Large_(131072)             1999           3119038 ns/op         672.37 MB/s     4194304 B/op          2 allocs/op
BenchmarkScientificFFT/Huge_(4194304)               22         243308959 ns/op         275.82 MB/s   134217728 B/op          2 allocs/op
BenchmarkScientificFFT/Massive_(16777216)            6        1015154183 ns/op         264.43 MB/s   536870912 B/op          2 allocs/op
BenchmarkFFT/Tiny_(4)                        780559590              7.72 ns/op        8286.02 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)                       4904836              1240 ns/op        1652.25 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)                       98176             62188 ns/op        1053.83 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)                       2062           3082183 ns/op         680.41 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)                         25         234863464 ns/op         285.74 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Massive_(16777216)                      6         988269633 ns/op         271.62 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Tiny_(4)                       400573915              14.8 ns/op        4334.42 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Small_(128)                      4230246              1433 ns/op        1429.33 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Medium_(4096)                      86562             68601 ns/op         955.32 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Large_(131072)                      1790           3320100 ns/op         631.65 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Huge_(4194304)                        22         257771536 ns/op         260.34 MB/s           0 B/op          0 allocs/op
BenchmarkIFFT/Massive_(16777216)                     5        1131598460 ns/op         237.22 MB/s           0 B/op          0 allocs/op
```

[travis-ci-badge]:   https://api.travis-ci.org/argusdusty/gofft.svg?branch=master
[travis-ci]:         https://api.travis-ci.org/argusdusty/gofft
[godoc-badge]:       https://godoc.org/github.com/argusdusty/gofft?status.svg
[godoc]:             https://godoc.org/github.com/argusdusty/gofft
[report-card-badge]: https://goreportcard.com/badge/github.com/argusdusty/gofft
[report-card]:       https://goreportcard.com/report/github.com/argusdusty/gofft