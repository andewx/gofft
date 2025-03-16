package fft

import (
	"math"
	"math/rand"
	"testing"
)

func TestIsPow2(t *testing.T) {
	// 1. Test all powers of 2 up to 2^63
	for i := 0; i < 64; i++ {
		x := 1 << uint64(i)
		r := IsPow2(x)
		if r != true {
			t.Errorf("IsPow2(%d), got: %t, expected: %t", x, r, true)
		}
	}

	// 2. Test all non-powers of 2 up to 2^15
	n := 1
	for x := 0; x < (1 << 16); x++ {
		if x == n {
			n <<= 1
			continue
		}
		r := IsPow2(x)
		if r != false {
			t.Errorf("IsPow2(%d), got: %t, expected: %t", x, r, false)
		}
	}
}

func TestNextPow2(t *testing.T) {
	// 0. Test n=0 returns 1
	r := NextPow2(0)
	if r != 1 {
		t.Errorf("NextPow2(0), got: %d, expected: 1", r)
	}
	for i := 0; i < 63; i++ {
		// 1. Test all powers of 2 up to 2^62
		x := 1 << uint32(i)
		r := NextPow2(x)
		if r != x {
			t.Errorf("NextPow2(%d), got: %d, expected: %d", x, r, x)
		}
		// 2. Test powers of 2 plus one
		r = NextPow2(x + 1)
		if r != 2*x {
			t.Errorf("NextPow2(%d+1), got: %d, expected: %d", x, r, 2*x)
		}
		// 3. Test random number between here and next power of 2
		if x > 1 {
			n := rand.Intn(x-1) + 1
			r = NextPow2(x + n)
			if r != 2*x {
				t.Errorf("NextPow2(%d+%d), got: %d, expected: %d", x, n, r, 2*x)
			}
		}
	}
}

func checkZeroPadding(t *testing.T, x1, x2 []complex128, N1, N2 int) {
	if len(x1) != N1 {
		t.Errorf("ZeroPad old array length, got: %d, expected: %d", len(x1), N1)
	}
	if len(x2) != N2 {
		t.Errorf("ZeroPad new array length, got: %d, expected: %d", len(x2), N2)
	}
	for j := 0; j < N1; j++ {
		if x1[j] != x2[j] {
			t.Errorf("ZeroPad copied section, got: x2[j] = %v, expected: x2[j] = %v", x2[j], x1[j])
		}
	}
	for j := N1; j < N2; j++ {
		if x2[j] != 0 {
			t.Errorf("ZeroPad padded section, got: x2[j] = %v, expected: x2[j] = %v", x2[j], 0)
		}
	}
}

func TestZeroPad(t *testing.T) {
	for i := 0; i < 100; i++ {
		// Test random lengths between 0 and 10000, and random paddings between 0 and 1000
		N1 := rand.Intn(10000)
		N2 := N1 + rand.Intn(1000)
		x1 := complexRand(N1)
		x2 := ZeroPad(x1, N2)
		checkZeroPadding(t, x1, x2, N1, N2)
	}
}

func TestZeroPadToNextPow2(t *testing.T) {
	// 0. Test n=0 returns [0]
	r := ZeroPadToNextPow2(nil)
	if len(r) != 1 {
		t.Errorf("len(ZeroPadToNextPow2(nil)), got: %d, expected: 1", len(r))
	}
	for i := 0; i < 17; i++ {
		// 1. Test powers of 2 up to 2^16
		N1 := 1 << uint32(i)
		x1 := complexRand(N1)
		x2 := ZeroPadToNextPow2(x1)
		checkZeroPadding(t, x1, x2, N1, N1)
		// 2. Test powers of 2 plus one
		x1 = complexRand(N1 + 1)
		x2 = ZeroPadToNextPow2(x1)
		checkZeroPadding(t, x1, x2, N1+1, 2*N1)
		// 3. Test random number between here and next power of 2
		if N1 > 1 {
			n := rand.Intn(N1-1) + 1
			x1 = complexRand(N1 + n)
			x2 = ZeroPadToNextPow2(x1)
			checkZeroPadding(t, x1, x2, N1+n, 2*N1)
		}
	}
}

func TestFloat64ToComplex128Array(t *testing.T) {
	// Test random arrays of length 0 to 1000
	for i := 0; i < 1000; i++ {
		a := floatRand(i)
		b := Float64ToComplex128Array(a)
		if len(a) != len(b) {
			t.Errorf("Float64ToComplex128Array, got: len(b) = %v, expected: len(b) = %v", len(b), len(a))
		}
		for j := 0; j < i; j++ {
			if a[j] != real(b[j]) {
				t.Errorf("Float64ToComplex128Array, got: real(b[j]) = %v, expected: real(b[j]) = %v", real(b[j]), a[j])
			}
			if imag(b[j]) != 0 {
				t.Errorf("Float64ToComplex128Array, got: imag(b[j]) = %v, expected: imag(b[j]) = 0", imag(b[j]))
			}
		}
	}
}

func TestComplex128ToFloat64Array(t *testing.T) {
	// Test random arrays of length 0 to 1000
	for i := 0; i < 1000; i++ {
		a := complexRand(i)
		b := Complex128ToFloat64Array(a)
		if len(a) != len(b) {
			t.Errorf("Complex128ToFloat64Array, got: len(b) = %v, expected: len(b) = %v", len(b), len(a))
		}
		for j := 0; j < i; j++ {
			if real(a[j]) != b[j] {
				t.Errorf("Complex128ToFloat64Array, got: b[j] = %v, expected: b[j] = %v", b[j], real(a[j]))
			}
		}
	}
}

func TestRoundFloat64Array(t *testing.T) {
	// Test random arrays of length 0 to 1000
	for i := 0; i < 1000; i++ {
		a := floatRand(i)
		b := make([]float64, i)
		copy(b, a)
		RoundFloat64Array(b)
		for j := 0; j < i; j++ {
			if math.Round(a[j]) != b[j] {
				t.Errorf("RoundFloat64Array, got: math.Round(a[j]) = %v, expected: math.Round(a[j]) = %v", math.Round(a[j]), b[j])
			}
		}
	}
}
