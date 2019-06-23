// Package gofft provides a fast discrete Fourier transformation algorithm.
//
// Implemented is the 1-dimensional DFT of complex input data
// for with input lengths which are powers of 2.
//
// The algorithm is non-recursive, works in-place overwriting
// the input array, and requires O(1) additional space.
package gofft

import (
	"fmt"
	"math"
	"sync"
)

var (
	prepareLock sync.RWMutex
	factorsMap  = map[int][]complex128{}
	permMap     = map[int][]int{}
)

// Prepare precomputes values used for FFT on a vector of length N.
// N must be a perfect power of 2, otherwise this will return an error.
func Prepare(N int) error {
	if !IsPow2(N) {
		return fmt.Errorf("Input dimension must be power of 2, is: %d", N)
	}
	prepareLock.RLock()
	if _, ok := factorsMap[N]; ok {
		prepareLock.RUnlock()
		// Already prepared, no need to do anything
		return nil
	}
	prepareLock.RUnlock()
	prepareLock.Lock()
	defer prepareLock.Unlock()
	factorsMap[N] = roots(N)
	permMap[N] = permutationIndex(N)
	return nil
}

// FFT implements the fast Fourier transform.
// This is done in-place (modifying the input array).
// Requires O(1) additional memory.
// len(x) must be a perfect power of 2, otherwise this will return an error.
func FFT(x []complex128) error {
	N, factors, perm, err := getVars(x)
	if err != nil {
		return err
	}
	fft(x, N, factors, perm)
	return nil
}

// IFFT implements the inverse fast Fourier transform.
// This is done in-place (modifying the input array).
// Requires O(1) additional memory.
// len(x) must be a perfect power of 2, otherwise this will return an error.
func IFFT(x []complex128) error {
	N, factors, perm, err := getVars(x)
	if err != nil {
		return err
	}
	ifft(x, N, factors, perm)
	return nil
}

// Pre-load the fft variables for later use.
func getVars(x []complex128) (N int, factors []complex128, perm []int, err error) {
	N = len(x)
	err = Prepare(N)
	factors = factorsMap[N]
	perm = permMap[N]
	return
}

// fft does the actual work for FFT
func fft(x []complex128, N int, factors []complex128, perm []int) {
	// Handle small N quickly
	switch N {
	case 1:
		return
	case 2:
		f := factors[0] * x[1]
		x[0], x[1] = x[0]+f, x[0]-f
		return
	case 4:
		x[1], x[2] = x[2], x[1]
		f := factors[0] * x[1]
		x[0], x[1] = x[0]+f, x[0]-f
		f = factors[0] * x[3]
		x[2], x[3] = x[2]+f, x[2]-f
		f = factors[0] * x[2]
		x[0], x[2] = x[0]+f, x[0]-f
		f = factors[1] * x[3]
		x[1], x[3] = x[1]+f, x[1]-f
		return
	}
	// Reorder the input array.
	permute(x, perm, N)
	// Butterfly
	s := N
	for n := 1; n < N; n <<= 1 {
		s >>= 1
		for o := 0; o < N; o += (n << 1) {
			for k := 0; k < n; k++ {
				i := k + o
				f := factors[k*s] * x[i+n]
				x[i], x[i+n] = x[i]+f, x[i]-f
			}
		}
	}
}

// ifft does the actual work for IFFT
func ifft(x []complex128, N int, factors []complex128, perm []int) {
	// Reverse the input vector
	for i := 1; i < N/2; i++ {
		j := N - i
		x[i], x[j] = x[j], x[i]
	}

	// Do the transform.
	fft(x, N, factors, perm)

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
	factors := make([]complex128, N/2)
	for n := 0; n < N/2; n++ {
		s, c := math.Sincos(-2.0 * math.Pi * float64(n) / float64(N))
		factors[n] = complex(c, s)
	}
	return factors
}
