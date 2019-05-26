package gofft

import (
	ktyefft "github.com/ktye/fft"
	"math/cmplx"
	"math/rand"
	"reflect"
	"testing"
)

// Run the benchmark against direct implementations with:
// go test -bench=.
//
// A benchmark against other more sophisticated implementations would be nice.

func complexRand(N int) []complex128 {
	x := make([]complex128, N)
	for i := 0; i < N; i++ {
		x[i] = complex(2.0*rand.Float64()-1.0, 2.0*rand.Float64()-1.0)
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

func BenchmarkSlow000(t *testing.B) {
	N := 1 << 13

	slow := slow{}
	x := complexRand(N)

	for i := 0; i < t.N; i++ {
		_ = slow.Transform(copyVector(x))
	}
}

func BenchmarkSlowPre(t *testing.B) {
	N := 1 << 13

	slowPre := newSlowPre(N)
	x := complexRand(N)

	for i := 0; i < t.N; i++ {
		_ = slowPre.Transform(copyVector(x))
	}
}

func BenchmarkKtye(t *testing.B) {
	N := 1 << 13

	f, err := ktyefft.New(N)
	if err != nil {
		t.Errorf("fft.New error: %v", err)
	}
	x := complexRand(N)

	for i := 0; i < t.N; i++ {
		f.Transform(copyVector(x))
	}
}

func BenchmarkFast(t *testing.B) {
	N := 1 << 13

	err := Prepare(N)
	if err != nil {
		t.Errorf("Prepare error: %v", err)
	}
	x := complexRand(N)

	for i := 0; i < t.N; i++ {
		err := FFT(copyVector(x))
		if err != nil {
			t.Errorf("FFT error: %v", err)
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

func TestIsPow2(t *testing.T) {
	// 1. Test all powers of 2 up to 2^30
	for i := 0; i < 31; i++ {
		x := 1 << uint32(i)
		r := isPow2(x)
		if r != true {
			t.Errorf("isPow2(%d), got: %t, expected: %t", x, r, true)
		}
	}

	// 2. Test all non-powers of 2 up to 2^15
	n := 1
	for x := 0; x < (1 << 15); x++ {
		if x == n {
			n <<= 1
			continue
		}
		r := isPow2(x)
		if r != false {
			t.Errorf("isPow2(%d), got: %t, expected: %t", x, r, false)
		}
	}
}

func TestBitReversal(t *testing.T) {
	tab := [][]int{
		[]int{0},
		[]int{0, 1},
		[]int{0, 2, 1, 3},
		[]int{0, 4, 2, 6, 1, 5, 3, 7},
		[]int{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15},
	}
	for i := 0; i < len(tab); i++ {
		got := permutationIndex(1 << uint32(i))
		expect := tab[i]
		if !reflect.DeepEqual(got, expect) {
			t.Errorf("%d expected: %v, got: %v\n", i, expect, got)
		}
	}
}
