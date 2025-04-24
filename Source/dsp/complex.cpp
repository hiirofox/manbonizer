#include "complex.h"

complex_f32_t cadd(complex_f32_t a, complex_f32_t b) { return complex_f32_t{ a.re + b.re, a.im + b.im }; }
complex_f32_t csub(complex_f32_t a, complex_f32_t b) { return complex_f32_t{ a.re - b.re, a.im - b.im }; }
complex_f32_t cmul(complex_f32_t a, complex_f32_t b) { return complex_f32_t{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re }; }
complex_f32_t cmulreal(complex_f32_t a, float b) { return complex_f32_t{ a.re * b, a.im * b }; }

complex_f32_t cdiv(complex_f32_t a, complex_f32_t b)
{
	float x = b.re * b.re + b.im * b.im;
	return complex_f32_t{ (a.re * b.re + a.im * b.im) / x,(a.im * b.re - a.re * b.im) / x };
}

complex_f64_t cadd64(complex_f64_t a, complex_f64_t b) { return complex_f64_t{ a.re + b.re, a.im + b.im }; }
complex_f64_t csub64(complex_f64_t a, complex_f64_t b) { return complex_f64_t{ a.re - b.re, a.im - b.im }; }
complex_f64_t cmul64(complex_f64_t a, complex_f64_t b) { return complex_f64_t{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re }; }
complex_f64_t cmulreal64(complex_f64_t a, double b) { return complex_f64_t{ a.re * b, a.im * b }; }