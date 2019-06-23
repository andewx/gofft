package gofft

import (
	"math"
)

// IsPow2 returns true if N is a perfect power of 2 (1, 2, 4, 8, ...)
// and false otherwise. Only works up to 2^30
// Algorithm from: https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
func IsPow2(N int) bool {
	if N == 0 {
		return false
	}
	return (N & (N - 1)) == 0
}

// NextPow2 returns the smallest power of 2 >= N. Only works up to 2^30
// Algorithm from: https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
func NextPow2(N int) int {
	if N == 0 {
		return 1
	}
	N -= 1
	N |= N >> 1
	N |= N >> 2
	N |= N >> 4
	N |= N >> 8
	N |= N >> 16
	return N + 1
}

// ZeroPad pads x with 0s at the end into a new array of length N.
// This does not alter x, and creates an entirely new array.
// This should only be used as a convience function, and isn't meant for performance.
// You should call this as few times as possible since it does potentially large allocations.
func ZeroPad(x []complex128, N int) []complex128 {
	y := make([]complex128, N)
	copy(y, x)
	return y
}

// ZeroPadToNextPow2 pads x with 0s at the end into a new array of length 2^N >= len(x)
// This does not alter x, and creates an entirely new array.
// This should only be used as a convience function, and isn't meant for performance.
// You should call this as few times as possible since it does potentially large allocations.
func ZeroPadToNextPow2(x []complex128) []complex128 {
	N := NextPow2(len(x))
	y := make([]complex128, N)
	copy(y, x)
	return y
}

// Float64ToComplex128Array converts a float64 array to the equivalent complex128 array
// using an imaginary part of 0.
func Float64ToComplex128Array(x []float64) []complex128 {
	y := make([]complex128, len(x))
	for i, v := range x {
		y[i] = complex(v, 0)
	}
	return y
}

// Complex128ToFloat64Array converts a complex128 array to the equivalent float64 array
// taking only the real part.
func Complex128ToFloat64Array(x []complex128) []float64 {
	y := make([]float64, len(x))
	for i, v := range x {
		y[i] = real(v)
	}
	return y
}

// RoundFloat64Array calls math.Round on each entry in x, changing the array in-place
func RoundFloat64Array(x []float64) {
	for i, v := range x {
		x[i] = math.Round(v)
	}
}
