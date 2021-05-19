#pragma once
#pragma once
// the basic vector in the 
#include <cmath>
using namespace std;
struct Vector3 {
	double x;
	double y;
	double z;
	// default method;
	Vector3() {
		x = double(0.0);
		y = double(0.0);
		z = double(0.0);
	}
	Vector3(double _x, double _y, double _z) {
		x = _x;
		y = _y;
		z = _z;
	}
	double dot(Vector3 vec1)const {
		return vec1.x * x + vec1.y * y + vec1.z * z;
	}
	Vector3& Normalize() {
		double distance = x * x + y * y + z * z;
		if (distance > 0) {
			double sqrtd = double(1.0) / sqrt(distance);
			x = x * sqrtd;
			y = y * sqrtd;
			z = z * sqrtd;
		}
		return *this;
	}

	Vector3& cross(Vector3 vec1, Vector3 vec2) {
		x = vec1.y * vec2.z - vec1.z * vec2.y;
		y = vec1.z * vec2.x - vec1.x * vec2.z;
		z = vec1.x * vec2.y - vec1.y * vec2.x;
		return *this;
	}
	//static Vector3 Cross(Vector3 vec1, Vector3 vec2) {
	//	return Vector3(vec1.y * vec2.z - vec1.z * vec2.y,
	//		vec1.z * vec2.x - vec1.x * vec2.z,
	//		vec1.x * vec2.y - vec1.y * vec2.x);
	//}

	Vector3 operator*(double value) {
		return Vector3(x *value, y * value, z * value);
	}


	bool operator==(Vector3 vec1) {
		return (x == vec1.x && y == vec1.y && z == vec1.z);
	}

	bool operator!=(Vector3 vec1) {
		return (x != vec1.x || y != vec1.y || z != vec1.z);
	}

	Vector3 operator-(Vector3 vec1) {
		return Vector3(x - vec1.x, y - vec1.y, z - vec1.z);
	}

	Vector3 operator+(Vector3 vec1) {
		return Vector3(x + vec1.x, y + vec1.y, z + vec1.z);
	}
	double len2() const { return (x*x + y * y + z * z); }
	double len() const { return std::sqrt(len2()); }
};