package gofft

import (
	"math"
	"math/cmplx"
	"math/rand"
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
	for i := 0; i < 64; i++ {
		x := complexRand(i)
		for j := 0; j < 64; j++ {
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
				if e := cmplx.Abs(r1[k] - r2[k]); e > 1e-9 {
					t.Errorf("slowConvolve and Convolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v", k, r1[k], k, r2[k], e)
				}
			}
		}
	}
}

func TestFastConvolve(t *testing.T) {
	// Test FastConvolve of zero inputs returns nil
	x := complexRand(0)
	y := complexRand(0)
	err := FastConvolve(x, y)
	if err != nil {
		t.Errorf("FastConvolve on empty inputs returned error: %v", err)
	}
	// Test FastConvolve of non-powers of 2 returns InputSizeError
	for N := 3; N < 101; N += 2 {
		x := complexRand(5)
		y := complexRand(5)
		err := FastConvolve(x, y)
		if err == nil {
			t.Errorf("FastConvolve on non-power of 2 input size %d didn't return error", N)
		}
		switch e := err.(type) {
		case *InputSizeError:
		default:
			t.Errorf("FastConvolve on non-power of 2 input size %d returned incorrect error type: %v", N, e)
		}
	}
	// Test FastConvolve of different vector sizes returns InputSizeError
	x = complexRand(4)
	y = complexRand(8)
	err = FastConvolve(x, y)
	if err == nil {
		t.Errorf("FastConvolve on differing input sizes didn't return error")
	}
	switch e := err.(type) {
	case *InputSizeError:
	default:
		t.Errorf("FastConvolve on differing input sizes returned incorrect error type: %v", e)
	}
	// Test FastConvolve(x, y) == slowConvolve(x, y)
	for i := 1; i < 128; i++ {
		N := NextPow2(2 * i)
		x := complexRand(i)
		x = ZeroPad(x, N)
		y := complexRand(i)
		y = ZeroPad(y, N)
		r1 := slowConvolve(x, y)
		err := FastConvolve(x, y)
		if err != nil {
			t.Error(err)
		}
		for j := 0; j < 2*i-1; j++ {
			if e := cmplx.Abs(r1[j] - x[j]); e > 1e-9 {
				t.Errorf("slowConvolve and FastConvolve differ: r1[%d]=%v, x[%d]=%v, diff=%v", j, r1[j], j, x[j], e)
			}
		}
		for j := 2*i - 1; j < N; j++ {
			if e := cmplx.Abs(x[j]); e > 1e-9 {
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
	// Test MultiConvolve() == nil
	x, err := MultiConvolve()
	if err != nil {
		t.Errorf("MultiConvolve() returned error: %v", err)
	}
	if len(x) != 0 {
		t.Errorf("MultiConvolve() returned non-empty result: %v", x)
	}
	// Test MultiConvolve(nil) == nil
	x, err = MultiConvolve(nil)
	if err != nil {
		t.Errorf("MultiConvolve(nil) returned error: %v", err)
	}
	if len(x) != 0 {
		t.Errorf("MultiConvolve(nil) returned non-empty result: %v", x)
	}
	// Test MultiConvolve(X...) == slowMultiConvolve(X)
	for i := 1; i < 25; i++ {
		X := make([][]complex128, i)
		for j := 1; j < 25; j++ {
			// Error propagates on the order of j^i
			errorThreshold := math.Pow(float64(j), float64(i)-1) * 1e-10
			// i arrays of length random(1, j)
			for k := 0; k < i; k++ {
				X[k] = complexRand(rand.Intn(j) + 1)
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
				if e := cmplx.Abs(r1[k] - r2[k]); e > errorThreshold {
					t.Errorf("slowMultiConvolve and MultiConvolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v, i=%d, j=%d", k, r1[k], k, r2[k], e, i, j)
				}
			}
		}
	}
}

func TestFastMultiConvolve(t *testing.T) {
	// Test non-power of 2 number of arrays
	err := FastMultiConvolve(make([]complex128, 5), 4, false)
	checkIsInputSizeError(t, "FastMultiConvolve(make([]complex128, 5), 4, false)", err)
	err = FastMultiConvolve(make([]complex128, 4), 3, false)
	checkIsInputSizeError(t, "FastMultiConvolve(make([]complex128, 4), 3, false)", err)
	err = FastMultiConvolve(make([]complex128, 4), 8, false)
	checkIsInputSizeError(t, "FastMultiConvolve(make([]complex128, 4), 8, false)", err)
	err = FastMultiConvolve(make([]complex128, 12), 4, false)
	checkIsInputSizeError(t, "FastMultiConvolve(make([]complex128, 12), 4, false)", err)
	// Test FastMultiConvolve(X) == slowMultiConvolve(X)
	for i := 1; i < 25; i++ {
		X1 := make([][]complex128, i)
		n := NextPow2(i)
		for j := 1; j < 25; j++ {
			// Error propagates on the order of j^i
			errorThreshold := math.Pow(float64(j), float64(i)-1) * 1e-10
			// i arrays of length j
			m := NextPow2(2 * j)
			X2 := make([]complex128, n*m)
			for k := 0; k < i; k++ {
				X1[k] = complexRand(j)
				copy(X2[m*k:m*(k+1)], X1[k])
			}
			for k := i; k < n; k++ {
				X2[m*k] = 1.0
			}
			X3 := make([]complex128, n*m)
			copy(X3, X2)
			r1 := slowMultiConvolve(X1)
			err := FastMultiConvolve(X2, m, false)
			if err != nil {
				t.Error(err)
			}
			err = FastMultiConvolve(X3, m, true)
			if err != nil {
				t.Error(err)
			}
			r2 := X2[:i*(j-1)+1]
			if len(r1) != len(r2) {
				t.Errorf("slowMultiConvolve and FastMultiConvolve differ in length: len(r1)=%d, len(r2)=%d", len(r1), len(r2))
			}
			r3 := X3[:i*(j-1)+1]
			if len(r2) != len(r3) {
				t.Errorf("FastMultiConvolve multithreaded difference in length: len(r2)=%d, len(r3)=%d", len(r2), len(r3))
			}
			for k := 0; k < len(r1); k++ {
				if e := cmplx.Abs(r1[k] - r2[k]); e > errorThreshold {
					t.Errorf("slowMultiConvolve and FastMultiConvolve differ: r1[%d]=%v, r2[%d]=%v, diff=%v, i=%d, j=%d", k, r1[k], k, r2[k], e, i, j)
				}
				if e := cmplx.Abs(r2[k] - r3[k]); e > errorThreshold {
					t.Errorf("FastMultiConvolve multithreaded difference: r2[%d]=%v, r3[%d]=%v, diff=%v, i=%d, j=%d", k, r2[k], k, r3[k], e, i, j)
				}
			}
		}
	}
}

func BenchmarkConvolve(b *testing.B) {
	for _, bm := range benchmarks {
		x := complexRand(bm.size)
		y := complexRand(bm.size)

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 32))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				Convolve(x, y)
			}
		})
	}
}

func BenchmarkFastConvolve(b *testing.B) {
	for _, bm := range benchmarks {
		x := complexRand(bm.size)
		y := complexRand(bm.size)

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 32))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				FastConvolve(x, y)
			}
		})
	}
}

var (
	multiConvolveBenchmarks = []struct {
		size   int
		number int
		name   string
	}{
		{4, 4, "Tiny (4, 4)"},
		{4096, 4, "Small (4096, 4)"},
		{131072, 4, "Medium-Tiny (131072, 4)"},
		{4096, 128, "Medium-Small (4096, 4096)"},
		{128, 4096, "Medium-Medium (128, 128)"},
		{4, 131072, "Medium-Large (4, 131072)"},
		{100000, 5, "Large-Tiny (100000, 5)"},
		{4000, 125, "Large-Small (4000, 125)"},
		{125, 4000, "Large-Medium (125, 4000)"},
		{5, 100000, "Large-Large (5, 100000)"},
	}
)

func BenchmarkMultiConvolve(b *testing.B) {
	for _, bm := range multiConvolveBenchmarks {
		x := make([][]complex128, bm.number)
		for i := 0; i < bm.number; i++ {
			x[i] = complexRand(bm.size)
		}

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * bm.number * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				MultiConvolve(x...)
			}
		})
	}
}

func BenchmarkFastMultiConvolve(b *testing.B) {
	for _, bm := range multiConvolveBenchmarks {
		x := make([]complex128, 2*NextPow2(bm.size)*NextPow2(bm.number))
		for i := 0; i < len(x); i += 2 * NextPow2(bm.size) {
			copy(x[i:], complexRand(NextPow2(bm.size)))
		}

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * bm.number * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				FastMultiConvolve(x, 2*NextPow2(bm.size), false)
			}
		})
	}
}

func BenchmarkFastMultiConvolveParallel(b *testing.B) {
	for _, bm := range multiConvolveBenchmarks {
		x := make([]complex128, 2*NextPow2(bm.size)*NextPow2(bm.number))
		for i := 0; i < len(x); i += 2 * NextPow2(bm.size) {
			copy(x[i:], complexRand(NextPow2(bm.size)))
		}

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * bm.number * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				err := FastMultiConvolve(x, 2*NextPow2(bm.size), true)
				if err != nil {
					b.Error(err)
				}
			}
		})
	}
}
