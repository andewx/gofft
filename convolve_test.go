package gofft

import (
	"math"
	"math/cmplx"
	"testing"
)

func slowConvolve(x, y []complex128) []complex128 {
	if len(x) == 0 && len(y) == 0 {
		return nil
	}
	r := make([]complex128, len(x)+len(y)-1)
	for i := 0; i < len(x); i++ {
		for j := 0; j < len(y); j++ {
			r[i+j] += x[i] * y[j]
		}
	}
	return r
}

func TestConvolve(t *testing.T) {
	for i := 0; i < 100; i++ {
		x := complexRand(i)
		for j := 0; j < 100; j++ {
			y := complexRand(j)
			r1 := slowConvolve(x, y)
			r2, err := Convolve(x, y)
			if err != nil {
				t.Error(err)
			}
			if len(r1) != len(r2) {
				t.Errorf("slowConvolve and Convolve differ in length: len(r1)=%d, len(r2)=%d", len(r1), len(r2))
			}
			for k := 0; k < i+j-1; k++ {
				if e := cmplx.Abs(r1[k] - r2[k]); e > 1E-9 {
					t.Errorf("slowConvolve and Convolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v", k, r1[k], k, r2[k], e)
				}
			}
		}
	}
}

func TestFastConvolve(t *testing.T) {
	for i := 1; i < 500; i++ {
		N := nextPow2(2 * i)
		err := Prepare(N)
		if err != nil {
			t.Error(err)
		}
		x := complexRand(i)
		x = ZeroPad(x, N)
		y := complexRand(i)
		y = ZeroPad(y, N)
		r1 := slowConvolve(x, y)
		err = FastConvolve(x, y)
		if err != nil {
			t.Error(err)
		}
		for j := 0; j < 2*i-1; j++ {
			if e := cmplx.Abs(r1[j] - x[j]); e > 1E-9 {
				t.Errorf("slowConvolve and FastConvolve differ: r1[%d]=%v, x[%d]=%v, diff=%v", j, r1[j], j, x[j], e)
			}
		}
		for j := 2*i - 1; j < N; j++ {
			if e := cmplx.Abs(x[j]); e > 1E-9 {
				t.Errorf("FastConvolve failed to zero-pad x: got x[%d]=%v, expected x[%d]=%v, diff=%v", j, x[j], j, 0, e)
			}
		}
		for j := 0; j < N; j++ {
			if y[j] != 0 {
				t.Errorf("FastConvolve failed to erase y: got y[%d]=%v, expected y[%d]=%v", j, y[j], j, 0)
			}
		}
	}
}

func slowMultiConvolve(X [][]complex128) []complex128 {
	m := []complex128{1.0}
	for _, x := range X {
		m = slowConvolve(m, x)
	}
	return m
}

func TestMultiConvolve(t *testing.T) {
	for i := 1; i < 25; i++ {
		X := make([][]complex128, i)
		for j := 1; j < 25; j++ {
			// Error propagates on the order of i^j
			error_threshold := math.Pow(float64(j), float64(i)-1) * 1E-11
			// i arrays of length j
			for k := 0; k < i; k++ {
				X[k] = complexRand(j)
			}
			r1 := slowMultiConvolve(X)
			r2, err := MultiConvolve(X...)
			if err != nil {
				t.Error(err)
			}
			if len(r1) != len(r2) {
				t.Errorf("slowMultiConvolve and MultiConvolve differ in length: len(r1)=%d, len(r2)=%d", len(r1), len(r2))
			}
			for k := 0; k < len(r1); k++ {
				if e := cmplx.Abs(r1[k] - r2[k]); e > error_threshold {
					t.Errorf("slowMultiConvolve and MultiConvolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v, i=%d, j=%d", k, r1[k], k, r2[k], e, i, j)
				}
			}
		}
	}
}

func TestFastMultiConvolve(t *testing.T) {
	for i := 1; i < 25; i++ {
		X1 := make([][]complex128, i)
		n := nextPow2(i)
		for j := 1; j < 25; j++ {
			// Error propagates on the order of i^j
			error_threshold := math.Pow(float64(j), float64(i)-1) * 1E-11
			// i arrays of length j
			m := nextPow2(2 * j)
			X2 := make([]complex128, n*m)
			for k := 0; k < i; k++ {
				X1[k] = complexRand(j)
				copy(X2[m*k:m*(k+1)], X1[k])
			}
			for k := i; k < n; k++ {
				X2[m*k] = 1.0
			}
			r1 := slowMultiConvolve(X1)
			err := FastMultiConvolve(X2, m, false)
			if err != nil {
				t.Error(err)
			}
			r2 := X2[:i*(j-1)+1]
			if len(r1) != len(r2) {
				t.Errorf("slowMultiConvolve and FastMultiConvolve differ in length: len(r1)=%d, len(r2)=%d", len(r1), len(r2))
			}
			for k := 0; k < len(r1); k++ {
				if e := cmplx.Abs(r1[k] - r2[k]); e > error_threshold {
					t.Errorf("slowMultiConvolve and FastMultiConvolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v, i=%d, j=%d", k, r1[k], k, r2[k], e, i, j)
				}
			}
		}
	}
}
