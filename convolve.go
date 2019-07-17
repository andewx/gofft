package gofft

import (
	"fmt"
	"runtime"
	"sync"
)

// Convolve computes the discrete convolution of x and y using FFT.
// Pads x and y to the next power of 2 from len(x)+len(y)-1
func Convolve(x, y []complex128) ([]complex128, error) {
	if len(x) == 0 && len(y) == 0 {
		return nil, nil
	}
	n := len(x) + len(y) - 1
	N := NextPow2(n)
	x = ZeroPad(x, N)
	y = ZeroPad(y, N)
	err := FastConvolve(x, y)
	return x[:n], err
}

// FastConvolve computes the discrete convolution of x and y using FFT
// and stores the result in x, while erasing y (setting it to 0s).
// Since this does no allocations, x and y are assumed to already be 0-padded
// for at least half their length.
func FastConvolve(x, y []complex128) error {
	if len(x) == 0 && len(y) == 0 {
		return nil
	}
	if len(x) != len(y) {
		return fmt.Errorf("x and y must have the same length, given: %d, %d", len(x), len(y))
	}
	convolve(x, y)
	return nil
}

// MultiConvolve computes the discrete convolution of many arrays using a
// hierarchical FFT algorithm that successfully builds up larger convolutions.
// This requires allocating up to 4*N extra memory for appropriate 0-padding
// where N=sum(len(x) for x in X).
// Takes O(N*log(N)^2) run time and O(N) additional space.
//
// This is much slower and takes many more allocations than FastMultiConvolve
// below, but has a smart planner that handles disproportionate array sizes
// very well. If all your arrays are the same length, FastMultiConvolve
// will be much faster.
func MultiConvolve(X ...[]complex128) ([]complex128, error) {
	arraysByLength := map[int][][]complex128{}
	mx := 1
	returnLength := 1
	for _, x := range X {
		// Pad out each array to the next power of two after twice the length
		// Doubling the length gives a buffer-zone for convolve to write to
		n := NextPow2(2 * len(x))
		arraysByLength[n] = append(arraysByLength[n], ZeroPad(x, n))
		if n > mx {
			mx = n
		}
		returnLength += len(x) - 1
	}
	if returnLength <= 0 {
		return nil, nil
	}
	// For each successive power of 2, convolve the entries in pairs up to the
	// next power of 2
	for i := 1; i <= mx; i *= 2 {
		arrays := arraysByLength[i]
		if len(arrays) > 0 {
			if len(arraysByLength) == 1 {
				return multiConvolveSingleLevel(arrays, returnLength)
			}
			for j := 0; j < len(arrays); j += 2 {
				if j+1 < len(arrays) {
					// For every pair, convolve to a single array
					convolve(arrays[j], arrays[j+1])
				}
				// Pad out to the next power of 2
				arraysByLength[2*i] = append(arraysByLength[2*i], ZeroPad(arrays[j], 2*i))
				if 2*i > mx {
					// Increase the max length as necessary
					// Shouldn't be possible to reach this, thanks to multiConvolveSingleLevel.
					mx = 2 * i
				}
			}
		}
		// Trigger the garbage collector
		arraysByLength[i] = nil
		delete(arraysByLength, i)
	}
	// Shouldn't be possible to reach this, but compiler needs it just in case.
	return arraysByLength[mx][0][:returnLength], nil
}

func multiConvolveSingleLevel(arrays [][]complex128, returnLength int) ([]complex128, error) {
	// If this is the final level, no need for further allocations,
	// just convolve together and return
	if len(arrays) == 2 {
		convolve(arrays[0], arrays[1])
		return arrays[0][:returnLength], nil
	}
	if len(arrays) == 1 {
		return arrays[0][:returnLength], nil
	}
	// If everything has ended up on the same length,
	// just use FastMultiConvolve and return
	N := len(arrays[0])
	n2 := NextPow2(len(arrays))
	data := make([]complex128, n2*N)
	for j, array := range arrays {
		copy(data[N*j:], array)
	}
	for j := len(arrays); j < n2; j++ {
		data[N*j] = 1.0
	}
	err := FastMultiConvolve(data, N, false)
	return data[:returnLength], err
}

// FastMultiConvolve computes the discrete convolution of many arrays using a
// hierarchical FFT algorithm, and stores the result in the first section of
// the input, writing 0s to the remainder of the input
// This does no allocations, so the arrays must first be 0-padded out to the
// next power of 2 from sum of the lengths of the longest two arrays.
// Additionally, the number of arrays must be a power of 2
// X is the concatenated array of arrays, of length N (n*m)
// n is the length of the 0-padded arrays.
// multithread tells the algorithm to use goroutines,
// which can slow things down for small N.
// Takes O(N*log(N)^2) run time and O(1) additional space.
func FastMultiConvolve(X []complex128, n int, multithread bool) error {
	N := len(X)
	if N%n != 0 {
		return fmt.Errorf("X must be array of arrays each of length n, instead have len(X) %d not divisible by n (%d)", N, n)
	}
	if !IsPow2(n) {
		return fmt.Errorf("X must be array of arrays each of a power of 2 length, instead have length %d not a power of 2", n)
	}
	if !IsPow2(N / n) {
		return fmt.Errorf("X must be array of arrays of a power of 2 length, instead have length %d not a power of 2", N/n)
	}
	for ; n != N; n <<= 1 {
		n2 := n << 1
		if multithread {
			var wg sync.WaitGroup
			NumCPU := runtime.NumCPU()
			for j := 0; j < NumCPU; j++ {
				wg.Add(1)
				go func(j int) {
					defer wg.Done()
					s := (j * (N / n2)) / NumCPU
					e := ((j + 1) * (N / n2)) / NumCPU
					for i := s; i < e; i++ {
						convolve(X[i*n2:i*n2+n], X[i*n2+n:i*n2+n2])
					}
				}(j)
			}
			wg.Wait()
		} else {
			for i := 0; i < N; i += n2 {
				convolve(X[i:i+n], X[i+n:i+n2])
			}
		}
	}
	return nil
}

// convolve does the actual work of convolutions.
func convolve(x, y []complex128) {
	fft(x)
	fft(y)
	for i := 0; i < len(x); i++ {
		x[i] *= y[i]
		y[i] = 0
	}
	ifft(x)
}
