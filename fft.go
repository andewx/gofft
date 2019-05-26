// Package fft provides a fast discrete Fourier transformation algorithm.
//
// Implemented is the 1-dimensional DFT of complex input data
// for with input lengths which are powers of 2.
//
// The algorithm is non-recursive, works in-place overwriting
// the input array, and requires O(1) additional space.
//
// Before doing the transform on acutal data, allocate
// an FFT object with t := fft.New(N) where N is the
// length of the input array.
// Then multiple calls to t.Transform(x) can be done with
// different input vectors having the same length.
package gofft

// ALGORITHM
// Example of the alogrithm with N=8 (P=3)
// Butterfly diagram:
//
//      1st stage p=0               2nd state p=1              3rd stage p=2
// IN +------------------+     +--------------------+     +-----------------------+
//                     overwrite                  overwrite                       output
// x0 -\/- x0 + E20 * x1 -> x0 -\  /- x0 + E40 * x2 -> x0 -\     /- x0 + E80 * x4 -> x0
// x1 -/\- x0 + E21 * x1 -> x1 -\\//- x1 + E41 * x3 -> x1 -\\   //- x1 + E81 * x5 -> x1
//                               /\                         \\ //
// x2 -\/- x0 + E22 * x1 -> x2 -//\\- x0 + E42 * x2 -> x2 -\ \\/ /- x2 + E82 * x6 -> x2
// x3 -/\- x0 + E23 * x1 -> x3 -/  \- x1 + E43 * x3 -> x3 -\\/\\//- x3 + E83 * x7 -> x3
//                                                          \\/\\
// x4 -\/- x0 + E24 * x1 -> x4 -\  /- x4 + E44 * x6 -> x4 -//\\/\\- x0 + E84 * x4 -> x4
// x5 -/\- x0 + E25 * x1 -> x5 -\\//- x5 + E45 * x7 -> x5 -/ /\\ \- x1 + E85 * x5 -> x5
//                               /\                         // \\
// x6 -\/- x0 + E26 * x1 -> x6 -//\\- x4 + E46 * x6 -> x6 -//   \\- x2 + E86 * x6 -> x6
// x7 -/\- x0 + E27 * x1 -> x7 -/  \- x5 + E47 * x7 -> x7 -/     \- x3 + E87 * x7 -> x7
//
// Enk are the N complex roots of 1 which were precomputed in E[0]..E[N-1].
// The stride s is N/n, and the index in E is k*s mod N,
// so   E21 of the first stage  is E[1*8/2 mod 8] = E[4]. These are +/- 1 alternating.
// and  E45 of the second stage is E[5*8/4 mod 8] = E[2]. These are 1,-i,-1,i and again.
// E8k  are all the roots (with stride=1) in increasing order: E[k].
//
// Before starting with the first stage, the input array must be
// permutated by the bit-inverted order.

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
	if !isPow2(N) {
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
		if !isPow2(N) {
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
	// Reorder the input array.
	permute(x, perm)

	s := N // Stride
	for n := 1; n < N; n <<= 1 {
		s >>= 1
		for o := 0; o < N; o += (n << 1) {
			for k := 0; k < n; k++ {
				// Do the butterfly with index i and j.
				i := k + o
				j := i + n
				x[i], x[j] = x[i]+E[k*s]*x[j], x[i]+E[s*(k+n)]*x[j]
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
	for i := range x {
		x[i] *= invN
	}
}

// permutationIndex builds the bit-inverted index vector,
// which is needed to permutate the input data.
func permutationIndex(N int) []int {
	index := make([]int, N)
	index[0] = 0 // Initial sequence for N=1
	// For every next power of two, the
	// sequence is multiplied by 2 inplace.
	// Then the result is also appended to the
	// end and increased by one.
	for n := 1; n < N; n <<= 1 {
		for i := 0; i < n; i++ {
			index[i] <<= 1
			index[i+n] = index[i] + 1
		}
	}
	return index
}

// permutate permutes the input vector according to the permutation vector.
// Uses an in-place algorithm that on FFT permutation vectors runs in O(N) time.
//
// Probably possible to speed this up with more precomputation (takes ~3/2*N rn)
func permute(x []complex128, perm []int) {
	n := len(x)
	for i := 0; i < (n - 1); i++ {
		ind := perm[i]
		for ind < i {
			ind = perm[ind]
		}
		x[i], x[ind] = x[ind], x[i]
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
