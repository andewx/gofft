# gofft
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
github.com\argusdusty\gofft>go test -benchmem -bench=. -cpu=1 -benchtime 15s
goos: windows
goarch: amd64
pkg: github.com/argusdusty/gofft
BenchmarkSlowFFT000/Tiny_(4)           100000000               239 ns/op         267.27 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFT000/Small_(128)           100000            322866 ns/op           6.34 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFT000/Medium_(4096)            100         317680973 ns/op           0.21 MB/s       65536 B/op          1 allocs/op
BenchmarkSlowFFTPre/Tiny_(4)           100000000               194 ns/op         328.75 MB/s          64 B/op          1 allocs/op
BenchmarkSlowFFTPre/Small_(128)           200000            161562 ns/op          12.68 MB/s        2048 B/op          1 allocs/op
BenchmarkSlowFFTPre/Medium_(4096)            100         165308052 ns/op           0.40 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Tiny_(4)              200000000               119 ns/op         537.11 MB/s          64 B/op          1 allocs/op
BenchmarkKtyeFFT/Small_(128)             5000000              3687 ns/op         555.44 MB/s        2048 B/op          1 allocs/op
BenchmarkKtyeFFT/Medium_(4096)            200000            162729 ns/op         402.73 MB/s       65536 B/op          1 allocs/op
BenchmarkKtyeFFT/Large_(131072)             2000           9276712 ns/op         226.07 MB/s     2097152 B/op          1 allocs/op
BenchmarkGoDSPFFT/Tiny_(4)              10000000              3114 ns/op          20.55 MB/s         520 B/op         13 allocs/op
BenchmarkGoDSPFFT/Small_(128)            2000000             10926 ns/op         187.43 MB/s        5590 B/op         18 allocs/op
BenchmarkGoDSPFFT/Medium_(4096)           100000            173007 ns/op         378.80 MB/s      164381 B/op         23 allocs/op
BenchmarkGoDSPFFT/Large_(131072)            3000           8108654 ns/op         258.63 MB/s     5243448 B/op         28 allocs/op
BenchmarkGoDSPFFT/Huge_(4194304)              30         688724583 ns/op          97.44 MB/s   167772809 B/op         33 allocs/op
BenchmarkFFT/Tiny_(4)                  500000000              39.5 ns/op        1619.09 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Small_(128)                20000000              1678 ns/op        1220.27 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Medium_(4096)                200000             96871 ns/op         676.53 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Large_(131072)                 3000           6669503 ns/op         314.44 MB/s           0 B/op          0 allocs/op
BenchmarkFFT/Huge_(4194304)                   20         887077165 ns/op          75.65 MB/s           0 B/op          0 allocs/op
```
