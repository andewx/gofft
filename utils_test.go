package gofft

import (
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
