package gofft

import (
	"math/rand"
	"testing"
)

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

func TestZeroPad(t *testing.T) {
	for i := 0; i < 1000; i++ {
		N1 := rand.Intn(10000)
		N2 := N1 + rand.Intn(1000)
		x1 := complexRand(N1)
		x2 := ZeroPad(x1, N2)
		if len(x1) != N1 {
			t.Errorf("ZeroPad old array length, got: %d, expected: %d", len(x1), N1)
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
}
