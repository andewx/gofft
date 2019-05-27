package gofft

import (
	ktyefft "github.com/ktye/fft"
	dspfft "github.com/mjibson/go-dsp/fft"
	"math"
	"math/cmplx"
	"math/rand"
	"reflect"
	"testing"
)

// Slow is the simplest and slowest FFT transform, for testing purposes
type slow struct {
}

func (s slow) Transform(x []complex128) []complex128 {
	N := len(x)
	y := make([]complex128, N)
	for k := 0; k < N; k++ {
		for n := 0; n < N; n++ {
			phi := -2.0 * math.Pi * float64(k*n) / float64(N)
			s, c := math.Sincos(phi)
			y[k] += x[n] * complex(c, s)
		}
	}
	return y
}

// SlowPre uses a precomputed roots table.
type slowPre struct {
	E []complex128
	N int
}

func newSlowPre(N int) slowPre {
	var s slowPre
	s.E = roots(N)
	s.N = N
	return s
}

func (s slowPre) Transform(x []complex128) []complex128 {
	if len(x) != len(s.E) {
		panic("SlowPre has been initialized with a long length, or not at all.")
	}
	y := make([]complex128, s.N)
	for k := 0; k < s.N; k++ {
		for n := 0; n < s.N; n++ {
			y[k] += x[n] * s.E[k*n%s.N]
		}
	}
	return y
}

func floatRand(N int) []float64 {
	x := make([]float64, N)
	for i := 0; i < N; i++ {
		x[i] = rand.NormFloat64()
	}
	return x
}

func complexRand(N int) []complex128 {
	x := make([]complex128, N)
	for i := 0; i < N; i++ {
		x[i] = complex(rand.NormFloat64(), rand.NormFloat64())
	}
	return x
}

func copyVector(v []complex128) []complex128 {
	y := make([]complex128, len(v))
	copy(y, v)
	return y
}

func TestFFT(t *testing.T) {
	N := 1 << 10
	x := complexRand(N)
	slow := slow{}
	slowPre := newSlowPre(N)
	err := Prepare(N)
	if err != nil {
		t.Errorf("Prepare error: %v", err)
	}

	y1 := slow.Transform(copyVector(x))
	y2 := slowPre.Transform(copyVector(x))
	y3 := copyVector(x)
	err = FFT(y3)
	if err != nil {
		t.Errorf("FFT error: %v", err)
	}
	for i := 0; i < N; i++ {
		if e := cmplx.Abs(y1[i] - y2[i]); e > 1E-9 {
			t.Errorf("slow and slowPre differ: i=%d diff=%v\n", i, e)
		}
		if e := cmplx.Abs(y1[i] - y3[i]); e > 1E-9 {
			t.Errorf("slow and fast differ: i=%d diff=%v\n", i, e)
		}
	}
}

func TestIFFT(t *testing.T) {
	N := 256
	x := complexRand(N)
	err := Prepare(N)
	if err != nil {
		t.Errorf("Prepare error: %v", err)
	}
	y := copyVector(x)
	err = FFT(y)
	if err != nil {
		t.Errorf("FFT error: %v", err)
	}
	err = IFFT(y)
	if err != nil {
		t.Errorf("IFFT error: %v", err)
	}
	for i := range x {
		if e := cmplx.Abs(x[i] - y[i]); e > 1E-9 {
			t.Errorf("inverse differs %d: %v %v\n", i, x[i], y[i])
		}
	}
}

func TestPermutationIndex(t *testing.T) {
	tab := [][]int{
		[]int{0},
		[]int{0, 1},
		[]int{0, 2, 2, 3},
		[]int{0, 4, 2, 6, 4, 5, 6, 7},
		[]int{0, 8, 4, 12, 4, 10, 6, 14, 8, 9, 10, 13, 12, 13, 14, 15},
	}
	for i := 0; i < len(tab); i++ {
		got := permutationIndex(1 << uint32(i))
		expect := tab[i]
		if !reflect.DeepEqual(got, expect) {
			t.Errorf("%d expected: %v, got: %v\n", i, expect, got)
		}
	}
}

var (
	benchmarks = []struct {
		size int
		name string
	}{
		{4, "Tiny (4)"},
		{128, "Small (128)"},
		{4096, "Medium (4096)"},
		{131072, "Large (131072)"},
		{4194304, "Huge (4194304)"},
	}
)

func BenchmarkSlowFFT000(b *testing.B) {
	for _, bm := range benchmarks {
		if bm.size > 10000 {
			// Don't run sizes too big for slow
			continue
		}
		b.Run(bm.name, func(b *testing.B) {
			slow := slow{}
			x := complexRand(bm.size)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = slow.Transform(x)
			}
		})
	}
}

func BenchmarkSlowFFTPre(b *testing.B) {
	for _, bm := range benchmarks {
		if bm.size > 10000 {
			// Don't run sizes too big for slow
			continue
		}
		b.Run(bm.name, func(b *testing.B) {
			slowPre := newSlowPre(bm.size)
			x := complexRand(bm.size)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				_ = slowPre.Transform(x)
			}
		})
	}
}

func BenchmarkKtyeFFT(b *testing.B) {
	for _, bm := range benchmarks {
		if bm.size > 1048576 {
			// Max size for ktye's fft
			continue
		}
		b.Run(bm.name, func(b *testing.B) {
			f, err := ktyefft.New(bm.size)
			if err != nil {
				b.Errorf("fft.New error: %v", err)
			}
			x := complexRand(bm.size)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				f.Transform(x)
			}
		})
	}
}

func BenchmarkDSPFFT(b *testing.B) {
	for _, bm := range benchmarks {
		b.Run(bm.name, func(b *testing.B) {
			x := complexRand(bm.size)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				dspfft.FFT(x)
			}
		})
	}
}

func BenchmarkFFT(b *testing.B) {
	for _, bm := range benchmarks {
		b.Run(bm.name, func(b *testing.B) {
			err := Prepare(bm.size)
			if err != nil {
				b.Errorf("Prepare error: %v", err)
			}
			x := complexRand(bm.size)

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				err := FFT(x)
				if err != nil {
					b.Errorf("FFT error: %v", err)
				}
			}
		})
	}
}
