package fft

import (
	"testing"
)

func TestInputSizeError(t *testing.T) {
	e := &InputSizeError{"asdf", "qwer", 5}
	expect := "Size of asdf must be qwer, is: 5"
	got := e.Error()
	if expect != got {
		t.Errorf("InputSizeError.Error(), expected %s, got %s", expect, got)
	}
}

func checkIsInputSizeError(t *testing.T, context string, err error) {
	if err == nil {
		t.Errorf("%s didn't return error", context)
	}
	switch e := err.(type) {
	case *InputSizeError:
	default:
		t.Errorf("%s returned incorrect error type: %v", context, e)
	}
}
