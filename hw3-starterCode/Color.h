#pragma once

#include <cmath>;
using namespace std;
struct Color {
	double r;
	double g;
	double b;

	Color() {
		r = 0.0;
		g = 0.0;
		b = 0.0;
	}
	Color(double _r, double _g, double _b) {
		r = _r;
		g = _g;
		b = _b;
	}
	
	Color& clamp() {
		r = max(r, 0);
		r = min(r, 1);
		g = max(g, 0);
		g = min(g, 1);
		b = max(b, 0);
		b = min(b, 1);

		return *this;
	}

	Color operator+(Color color1) {
		r += color1.r;
		g += color1.g;
		b += color1.b;
		return *this;
	}

	Color& operator+=(Color color1) {
		r += color1.r;
		g += color1.g;
		b += color1.b;
		return *this;
	}

	Color operator*(double scaler) {
		return (Color(r * scaler, g * scaler, b* scaler));
	}

	Color operator/(double scaler) {
		if (scaler == 0) return *this;
		return (Color(r / scaler, g / scaler, b / scaler));
	}
};