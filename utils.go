package gofft

import (
	"math"
	"math/bits"
)

// IsPow2 returns true if N is a perfect power of 2 (1, 2, 4, 8, ...) and false otherwise.
// Algorithm from: https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
func IsPow2(N int) bool {
	if N == 0 {
		return false
	}
	return (uint64(N) & uint64(N-1)) == 0
}

// NextPow2 returns the smallest power of 2 >= N.
func NextPow2(N int) int {
	if N == 0 {
		return 1
	}
	return 1 << uint64(bits.Len64(uint64(N-1)))
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

// Complex64ToComplex128 converts a slice of complex64 to complex128
func Complex64ToComplex128(data []complex64) []complex128 {
	result := make([]complex128, len(data))
	for i := 0; i < len(data); i++ {
		result[i] = complex(float64(real(data[i])), float64(imag(data[i])))
	}
	return result
}

// Complex128ToComplex64 converts a slice of complex128 to complex64
func Complex128ToComplex64(data []complex128) []complex64 {
	result := make([]complex64, len(data))
	for i := 0; i < len(data); i++ {
		result[i] = complex(float32(real(data[i])), float32(imag(data[i])))
	}
	return result
}
