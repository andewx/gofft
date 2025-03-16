package fft

import (
	"math"
	"math/cmplx"
)

type Window int

const (
	Rectangular Window = iota
	Hanning
	Hamming
	Blackman
)

// ApplyWindow applies the specified window function to the input data
func ApplyWindow(x []complex128, window Window) []complex128 {
	n := len(x)

	for i := 0; i < n; i++ {
		var w float64
		switch window {
		case Rectangular:
			w = 1.0
		case Hanning:
			w = 0.5 * (1 - math.Cos(2*math.Pi*float64(i)/float64(n-1)))
		case Hamming:
			w = 0.54 - 0.46*math.Cos(2*math.Pi*float64(i)/float64(n-1))
		case Blackman:
			w = 0.42 - 0.5*math.Cos(2*math.Pi*float64(i)/float64(n-1)) +
				0.08*math.Cos(4*math.Pi*float64(i)/float64(n-1))
		}
		x[i] = complex(real(x[i])*w, imag(x[i])*w)
	}

	return x
}

// ApplyWindow64 applies the specified window function to the input data
func ApplyWindow64(x []complex64, window Window) []complex64 {
	n := len(x)

	for i := 0; i < n; i++ {
		var w float64
		switch window {
		case Rectangular:
			w = 1.0
		case Hanning:
			w = 0.5 * (1 - math.Cos(2*math.Pi*float64(i)/float64(n-1)))
		case Hamming:
			w = 0.54 - 0.46*math.Cos(2*math.Pi*float64(i)/float64(n-1))
		case Blackman:
			w = 0.42 - 0.5*math.Cos(2*math.Pi*float64(i)/float64(n-1)) +
				0.08*math.Cos(4*math.Pi*float64(i)/float64(n-1))
		}
		x[i] = complex(real(x[i])*float32(w), imag(x[i])*float32(w))
	}

	return x
}

// PowerSpectrumPrecision computes the power spectrum of the FFT result
func PowerSpectrumPrecision(x []complex128) []float64 {
	n := len(x)
	result := make([]float64, n)

	for i := 0; i < n; i++ {
		result[i] = cmplx.Abs(x[i]) * cmplx.Abs(x[i])
	}

	return result
}

// PowerSpectrum computes the power spectrum of the FFT result
func PowerSpectrum(x []complex64) []float32 {
	n := len(x)
	result := make([]float32, n)

	for i := 0; i < n; i++ {
		result[i] = real(x[i])*real(x[i]) + imag(x[i])*imag(x[i])
	}

	return result
}
