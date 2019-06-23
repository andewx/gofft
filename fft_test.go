package gofft

import (
	ktyefft "github.com/ktye/fft"
	dspfft "github.com/mjibson/go-dsp/fft"
	"math"
	"math/cmplx"
	"math/rand"
	"reflect"
	"runtime"
	"sync/atomic"
	"testing"
)

// Slow is the simplest and slowest FFT transform, for testing purposes
func slowFFT(x []complex128) []complex128 {
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
	err := Prepare(N)
	if err != nil {
		t.Errorf("Prepare error: %v", err)
	}

	y1 := slowFFT(copyVector(x))
	y2 := copyVector(x)
	err = FFT(y2)
	if err != nil {
		t.Errorf("FFT error: %v", err)
	}
	for i := 0; i < N; i++ {
		if e := cmplx.Abs(y1[i] - y2[i]); e > 1E-9 {
			t.Errorf("slowFFT and FFT differ: i=%d diff=%v\n", i, e)
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

func BenchmarkSlowFFT(b *testing.B) {
	for _, bm := range benchmarks {
		if bm.size > 10000 {
			// Don't run sizes too big for slow
			continue
		}
		x := complexRand(bm.size)

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				slowFFT(x)
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
		f, err := ktyefft.New(bm.size)
		if err != nil {
			b.Errorf("fft.New error: %v", err)
		}
		x := complexRand(bm.size)

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				f.Transform(x)
			}
		})
	}
}

func BenchmarkGoDSPFFT(b *testing.B) {
	for _, bm := range benchmarks {
		dspfft.EnsureRadix2Factors(bm.size)
		x := complexRand(bm.size)

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				dspfft.FFT(x)
			}
		})
	}
}

func BenchmarkFFT(b *testing.B) {
	for _, bm := range benchmarks {
		x := complexRand(bm.size)
		err := Prepare(bm.size)
		if err != nil {
			b.Errorf("Prepare error: %v", err)
		}

		b.Run(bm.name, func(b *testing.B) {
			b.SetBytes(int64(bm.size * 16))
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				FFT(x)
				if err != nil {
					b.Errorf("FFT error: %v", err)
				}
			}
		})
	}
}

func BenchmarkFFTParallel(b *testing.B) {
	for _, bm := range benchmarks {
		procs := runtime.GOMAXPROCS(0)
		x := complexRand(bm.size * procs)
		err := Prepare(bm.size)
		if err != nil {
			b.Errorf("Prepare error: %v", err)
		}

		b.Run(bm.name, func(b *testing.B) {
			var idx uint64
			b.SetBytes(int64(bm.size * 16))
			b.ResetTimer()
			b.RunParallel(func(pb *testing.PB) {
				i := int(atomic.AddUint64(&idx, 1) - 1)
				y := x[i*bm.size : (i+1)*bm.size]
				for pb.Next() {
					FFT(y)
				}
			})
		})
	}
}
