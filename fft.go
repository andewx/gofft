// Package fft provides a fast discrete Fourier transformation algorithm.
//
// Implemented is the 1-dimensional DFT of complex input data
// for with input lengths which are powers of 2.
//
// The algorithm is non-recursive, works in-place overwriting
// the input array, and requires O(1) additional space.
//
// Before doing the transform on acutal data, prepare the fft with
// t := fft.Prepare(N) where N is the length of the input array.
// Then multiple calls to fft.FFT(x) can be done with
// different input vectors having the same length.
package gofft

import (
	"fmt"
	"math"
)

var (
	Es    = map[int][]complex128{}
	perms = map[int][]int{}
)

// Prepare precomputes values used for FFT on a vector of length N.
// N must be a perfect power of 2, otherwise this will return an error.
func Prepare(N int) error {
	if !IsPow2(N) {
		return fmt.Errorf("Input dimension must be power of 2, is: %d", N)
	}
	if _, ok := Es[N]; ok {
		// Already prepared, no need to do anything
		return nil
	}
	Es[N] = roots(N)
	perms[N] = permutationIndex(N)
	return nil
}

// checkN tests N as being a valid FFT vector length.
// Returns an error if it isn't.
func checkN(N int) error {
	if _, ok := Es[N]; !ok {
		if !IsPow2(N) {
			return fmt.Errorf("Input dimension must be power of 2, is: %d", N)
		}
		return fmt.Errorf("FFT is not initialized for input dimension: %d, must initialize with Prepare(N) first", N)
	}
	return nil
}

// FFT implements the fast Fourier transform.
// This is done in-place (modifying the input array).
// Requires O(1) additional memory.
// len(x) must be a perfect power of 2, otherwise this will return an error.
// You must call Prepare(len(x)) before this, otherwise this will return an error.
func FFT(x []complex128) error {
	N, E, perm, err := getVars(x)
	if err != nil {
		return err
	}
	fft(x, N, E, perm)
	return nil
}

// IFFT implements the inverse fast Fourier transform.
// This is done in-place (modifying the input array).
// Requires O(1) additional memory.
// len(x) must be a perfect power of 2, otherwise this will return an error.
// You must call Prepare(len(x)) before this, otherwise this will return an error.
func IFFT(x []complex128) error {
	N, E, perm, err := getVars(x)
	if err != nil {
		return err
	}
	ifft(x, N, E, perm)
	return nil
}

// Pre-load the fft variables for later use.
func getVars(x []complex128) (N int, E []complex128, perm []int, err error) {
	N = len(x)
	E = Es[N]
	perm = perms[N]
	err = checkN(N)
	return
}

// fft does the actual work for FFT
func fft(x []complex128, N int, E []complex128, perm []int) {
	// Handle small N quickly
	switch N {
	case 1:
		return
	case 2:
		x[0], x[1] = x[0]+E[0]*x[1], x[0]+E[1]*x[1]
		return
	case 4:
		x[1], x[2] = x[2], x[1]
		x[0], x[1] = x[0]+E[0]*x[1], x[0]+E[2]*x[1]
		x[2], x[3] = x[2]+E[0]*x[3], x[2]+E[2]*x[3]
		x[0], x[2] = x[0]+E[0]*x[2], x[0]+E[2]*x[2]
		x[1], x[3] = x[1]+E[1]*x[3], x[1]+E[3]*x[3]
		return
	}
	// Reorder the input array.
	permute(x, perm, N)
	// Butterfly
	s := N
	E2 := E[N>>1:]
	for n := 1; n < N; n <<= 1 {
		s >>= 1
		for o := 0; o < N; o += (n << 1) {
			for k := 0; k < n; k++ {
				i := k + o
				x[i], x[i+n] = x[i]+E[k*s]*x[i+n], x[i]+E2[k*s]*x[i+n]
			}
		}
	}
}

// ifft does the actual work for IFFT
func ifft(x []complex128, N int, E []complex128, perm []int) {
	// Reverse the input vector
	for i := 1; i < N/2; i++ {
		j := N - i
		x[i], x[j] = x[j], x[i]
	}

	// Do the transform.
	fft(x, N, E, perm)

	// Scale the output by 1/N
	invN := complex(1.0/float64(N), 0)
	for i := 0; i < N; i++ {
		x[i] *= invN
	}
}

// permutationIndex builds the bit-inverted index vector,
// which is needed to permutate the input data.
func permutationIndex(N int) []int {
	index := make([]int, N)
	index[0] = 0 // Initial sequence for N=1
	// For every next power of two, the sequence is multiplied by 2 in-place.
	// Then the result is also appended to the end and increased by one.
	for n := 1; n < N; n <<= 1 {
		for i := 0; i < n; i++ {
			index[i] <<= 1
			index[i+n] = index[i] + 1
		}
	}
	// Re-arrange the permutation to just the necessary swaps
	for i := 1; i < N-1; i++ {
		ind := index[i]
		for ind < i {
			ind = index[ind]
		}
		index[i] = ind
	}
	return index
}

// permutate permutes the input vector according to the permutation vector.
// Uses an in-place algorithm that on FFT permutation vectors runs in O(N) time.
func permute(x []complex128, perm []int, N int) {
	// perm[0] is always 0, and perm[N-1] is always N-1, so skip those
	for i := 1; i < N-1; i++ {
		x[i], x[perm[i]] = x[perm[i]], x[i]
	}
}

// roots computes the complex-roots-of 1 table of length N.
func roots(N int) []complex128 {
	E := make([]complex128, N)
	for n := 0; n < N; n++ {
		s, c := math.Sincos(-2.0 * math.Pi * float64(n) / float64(N))
		E[n] = complex(c, s)
	}
	return E
}
