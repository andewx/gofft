package gofft

import (
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
					t.Errorf("slowConvolve and Convolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v\n", k, r1[k], k, r2[k], e)
				}
			}
		}
	}
}
