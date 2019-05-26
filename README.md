# gofft
radix-2 fast Fourier transform

Package gofft provides an fast discrete Fourier transformation algorithm in pure Go.

Implemented is the 1-dimensional DFT of complex input data for input lengths which are powers of 2.

The algorithm is non-recursive and works in-place overwriting the input array.

Requires constant (O(1)) additional space.

## Why?
Most existing FFT libraries in Go allocate temporary arrays with O(N) additional space. This is less-than-ideal when you have arrays of length of 2<sup>25</sup> or more, where you quickly end up allocating gigabytes of data and dragging down the FFT calculation to a halt.

I took [an existing](https://github.com/ktye/fft) FFT implementation in Go, cleaned and improved the code API and performance, and replaced the permutation step with an algorithm that works with no temp array. Performance was nearly doubled (see benchmarks in fft_test), is much simpler to use, and handles errors better.

Also included a new utility functions: `Pad(x, N)`, `FastConvolve(x, y)`, `Float64ToComplex128Array`, and `Complex128ToFloat64Array`.

## How?
```go
package main

import (
	"fmt"
	"github.com/argusdusty/gofft"
)

func main() {
	testArray := gofft.Float64ToComplex128Array([]float64{1, 2, 3, 4, 5, 6, 7, 8})
	gofft.Prepare(len(testArray))
	gofft.FFT(testArray)
	gofft.IFFT(testArray)
	fmt.Println(gofft.Complex128ToFloat64Array(testArray))
}
```
Outputs: `[1.0000000000000004 2.0000000000000004 3.0000000000000004 4 5 6 7 8]`