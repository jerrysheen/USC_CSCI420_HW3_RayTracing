/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: TaoShen
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include "Color.h"
#include <imageIO.h>
#include "Vector3.h"
#include <cmath>
 
#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
#define PI 3.1415926
//the field of view of the camera
#define fov 60.0
#define SMALL_VALUE 0.0000001
unsigned char buffer[HEIGHT][WIDTH][3];
bool enableAntialias = false;
bool enableReflection = false;
bool enableSoftShadow = false;
struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
int totalRefelctTime = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
Color phongShadingSphere(Sphere sphere, Light light, Vector3 intersectionPoint);
//
struct Ray{
	Vector3 origin;
	Vector3 direction;

	// constrctor.
	Ray() {}
	Ray(Vector3 o, Vector3 d) {
		origin = o;
		direction = d;
	}

	//
	bool isIntersectionSphere(const Sphere& s, Vector3& intersectionPoint) {

		double a = direction.dot(direction);
		//printf("%f   ", a);
		// xo_xc = xo
		Vector3 xo_xc = Vector3(origin.x - s.position[0], origin.y - s.position[1], origin.z - s.position[2]); 
		//printf("  %f, %f, %f,", s.position[0], s.position[1], s.position[2]);
		//printf("  %f, %f, %f,", direction.x, direction.y, direction.z);
		double b = 2.0 * direction.dot(xo_xc);
		double r = (double)s.radius;
		double c = xo_xc.dot(xo_xc) - r * r;

		//// b^2 - 4ac, a = 1;
		double discriminant = b * b - 4.0 *a * c;
		/*double t0 = 0.0, t1 = 0.0, t = 0.0;*/
		double t0, t1, t;
		t0 = t1 = t = 0.0;
		//printf("  %d, %d, %d,", xo_xc.x, xo_xc.y, xo_xc.z);
		//printf("  %f\n", discriminant);
		if (discriminant < -SMALL_VALUE) return false;
		if (discriminant < SMALL_VALUE) {
			discriminant = abs(discriminant);
			/*t0 = (-b) * 0.5;
			t1 = t0;*/
			t0 = (-b + std::sqrt(discriminant)) * 0.5 / a;; // (-b (+/-) sqrt(absD))/2a, but a = 1
			t1 = (-b - std::sqrt(discriminant)) * 0.5 / a;;
			t = t0 / a;
			//printf("  %f\n", discriminant);
		}
		else {
			//here I meet Loss of significance problem. and I find others follow: https://en.wikipedia.org/wiki/Loss_of_significance
			if (b >= 0) {
				t0 = (-b - sqrt(discriminant)) * (double)0.5;
				t1 = c / t0;
			}
			else {
				t0 = (-b + sqrt(discriminant)) * (double)0.5;
				t1 = c / t0;
			}
		}

		// assign t0 < t1;
		t = t1;
		if (t0 > t1) {
			t1 = t0;
			t0 = t;
		}

		if (t1 < 0) return false;
		if (t0 > 0) t = t0;
		else t = t1;
		intersectionPoint = origin + direction * t;
		return true;
	}

	bool isIntersectionTriangle(const Triangle& triangle, Vector3& intersectionPoint) {
		Vector3 vertex0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
		Vector3 vertex1(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
		Vector3 vertex2(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
		Vector3 normal;
		// the plane normal;
		normal.cross(vertex1 - vertex0, vertex2 - vertex0).Normalize();
		// if plane normal is penpendicular to the direction, or in a small value, 
		// then it means no intersection with the plane.
		// normalPlaneAngle = [a, b , c];
		double normalPlaneAngle = normal.dot(direction);
		if (abs(normalPlaneAngle) < SMALL_VALUE) return false;
		// ax + by + cx + d = 0
		// d = -(ax + by + cx) = - normal.dot(vertex 0)
		double d = -normal.dot(vertex0);
		// get t using the form s = (n.dot origin point + d) / n * direction
		double t = -(normal.dot(origin) + d) / normal.dot(direction);
		// t <= 0 means intersection point is behind or at the origin point.
		if (t <= 0) return false;
		intersectionPoint = origin + direction * t;
		// to find a point inside or out side the triangle check this:
		// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
		Vector3 a = vertex0 - intersectionPoint;
		Vector3 b = vertex1 - intersectionPoint;
		Vector3 c = vertex2 - intersectionPoint;
		Vector3 u, v, w;
		u.cross(b, c);
		v.cross(c, a);
		w.cross(a, b);

		if (u.dot(v) < 0.0f || u.dot(w) < 0.0f) return false;
		//printf("intersection!!");
		return true;
	}

};


Ray* rayList(double pixelX, double pixelY) {
	if (!enableAntialias) {
		double aspectRatio = (double)WIDTH / (double)HEIGHT;
		// define fov = 60;
		double tanFov = tan((60.0 * (double)PI / 180.0) * 0.5);
		Ray* ray = new Ray[1];

		// get the ratio:   [0, 1]
		double ratioX = (pixelX + 0.5) / (double)WIDTH;
		double ratioY = (pixelY + 0.5) / (double)HEIGHT;

		// in the NDC space its from [-1, 1]
		double NDCx = ratioX * 2 - 1;
		double NDCy = ratioY * 2 - 1;

		// transfrom it to fit the screen ratio 
		double x = NDCx * aspectRatio * tanFov;
		double y = NDCy * tanFov;
		/*double x = NDCx * imageAspectRatio * tanhalfFOV;
		double y = NDCy * tanhalfFOV;*/
		ray[0] = Ray(Vector3(0, 0, 0), Vector3(x, y, -1.0).Normalize());
		return ray;
	}
	else {
		// we create 4 rays and then average it color
		// define fov = 60;
		double aspectRatio = (double)WIDTH / (double)HEIGHT;
		double tanFov = tan((60.0 * (double)PI / 180.0) * 0.5);
		Ray* ray = new Ray[4];

		// get the ratio:   [0, 1]
		double ratioX = (pixelX + 0.25) / (double)WIDTH;
		double ratioY = (pixelY + 0.25) / (double)HEIGHT;
		// in the NDC space its from [-1, 1]
		double NDCx = ratioX * 2 - 1;
		double NDCy = ratioY * 2 - 1;
		// transfrom it to fit the screen ratio 
		double x = NDCx * aspectRatio * tanFov;
		double y = NDCy * tanFov;
		/*double x = NDCx * imageAspectRatio * tanhalfFOV;
		double y = NDCy * tanhalfFOV;*/
		ray[0] = Ray(Vector3(0, 0, 0), Vector3(x, y, -1.0).Normalize());

		// get the ratio:   [0, 1]
		ratioX = (pixelX + 0.25) / (double)WIDTH;
		ratioY = (pixelY + 0.75) / (double)HEIGHT;
		// in the NDC space its from [-1, 1]
		NDCx = ratioX * 2 - 1;
		NDCy = ratioY * 2 - 1;
		// transfrom it to fit the screen ratio 
		x = NDCx * aspectRatio * tanFov;
		y = NDCy * tanFov;
		/*double x = NDCx * imageAspectRatio * tanhalfFOV;
		double y = NDCy * tanhalfFOV;*/
		ray[1] = Ray(Vector3(0, 0, 0), Vector3(x, y, -1.0).Normalize());

		// get the ratio:   [0, 1]
		ratioX = (pixelX + 0.75) / (double)WIDTH;
		ratioY = (pixelY + 0.25) / (double)HEIGHT;
		// in the NDC space its from [-1, 1]
		NDCx = ratioX * 2 - 1;
		NDCy = ratioY * 2 - 1;
		// transfrom it to fit the screen ratio 
		x = NDCx * aspectRatio * tanFov;
		y = NDCy * tanFov;
		/*double x = NDCx * imageAspectRatio * tanhalfFOV;
		double y = NDCy * tanhalfFOV;*/
		ray[2] = Ray(Vector3(0, 0, 0), Vector3(x, y, -1.0).Normalize());

		// get the ratio:   [0, 1]
		ratioX = (pixelX + 0.75) / (double)WIDTH;
		ratioY = (pixelY + 0.75) / (double)HEIGHT;
		// in the NDC space its from [-1, 1]
		NDCx = ratioX * 2 - 1;
		NDCy = ratioY * 2 - 1;
		// transfrom it to fit the screen ratio 
		x = NDCx * aspectRatio * tanFov;
		y = NDCy * tanFov;
		/*double x = NDCx * imageAspectRatio * tanhalfFOV;
		double y = NDCy * tanhalfFOV;*/
		ray[3] = Ray(Vector3(0, 0, 0), Vector3(x, y, -1.0).Normalize());

		return ray;
	}

}

void clamp(double& value) {
	if (value > 1.0) value = 1.0;
	if (value < 0.0) value = 0.0;
}

// generate phong shading color by its normal and light color. 
Color phongShadingSphere(Sphere sphere, Light light, Vector3 intersectionPoint){
	//printf("%f, %f, %f", intersectionPoint.x, intersectionPoint.y, intersectionPoint.z);
	Color c;
	Vector3 sphereCenter = Vector3(sphere.position[0], sphere.position[1], sphere.position[2]);
	Vector3 lightSource = Vector3(light.position[0], light.position[1], light.position[2]);
	Vector3 N = Vector3(intersectionPoint - sphereCenter).Normalize();
	
	Vector3 L = Vector3(lightSource - intersectionPoint).Normalize();

	double LdotN = L.dot(N);
	clamp(LdotN);
	
	// if L dot N < 0 which means light and normal direction is opposite, so color 
	// should be zreo;

	// only can do vector3 * double instead of double * vector3;
	Vector3 R = (N * 2.0 * LdotN - L).Normalize();
	Vector3 V = (Vector3(0,0,0) - intersectionPoint).Normalize();
	double RdotV = R.dot(V);
	clamp(RdotV);
	Color kd(sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
	Color ks(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
	
	c.r = light.color[0] * (kd.r * LdotN) + ks.r * pow(RdotV, sphere.shininess);
	c.g = light.color[1] * (kd.g * LdotN) + ks.g * pow(RdotV, sphere.shininess);
	c.b = light.color[2] * (kd.b * LdotN) + ks.b * pow(RdotV, sphere.shininess);

	return c;
}

// generate phong shading color by its normal and light color. 
Color phongShadingTriangle(Triangle triangle, Light light, Vector3 intersectionPoint) {
	Vector3 c0 (triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
	Vector3 c1 (triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
	Vector3 c2 (triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
	Vector3 c = intersectionPoint;
	// interpolate normal with  barycentric coordinates method
	Vector3 temp0;
	double areaC0C1C2 = temp0.cross(c2 - c0, c1 - c0).len();
	// a = area(c c1 c2) / area(c0 c1 c2)
	double a = temp0.cross(c2 - c, c1 - c).len() / areaC0C1C2;
	double b = temp0.cross(c2 - c, c0 - c).len() / areaC0C1C2;
	//double y = temp0.cross(c1 - c, c0 - c).len() / areaC0C1C2;
	double y = 1 - a - b;

	// interpolate normal by p = ap1 + bp2 + yp3;
	Vector3 normal = Vector3(
		(a * triangle.v[0].normal[0]) + (b * triangle.v[1].normal[0]) + (y * triangle.v[2].normal[0]),
		(a * triangle.v[0].normal[1]) + (b * triangle.v[1].normal[1]) + (y * triangle.v[2].normal[1]),
		(a * triangle.v[0].normal[2]) + (b * triangle.v[1].normal[2]) + (y * triangle.v[2].normal[2])).Normalize();
	
	Vector3 lightSource(light.position[0], light.position[1], light.position[2]);
	Vector3 L = Vector3(lightSource - intersectionPoint).Normalize();

	double LdotN = L.dot(normal);
	clamp(LdotN);
	// if L dot N < 0 which means light and normal direction is opposite, so color 
	// should be zreo;
	// only can do vector3 * double instead of double * vector3;
	Vector3 R = (normal * 2.0 * LdotN - L).Normalize();
	// V = Eye point
	Vector3 E = (Vector3(0, 0, 0) - intersectionPoint).Normalize();
	double RdotE = R.dot(E);
	clamp(RdotE);

	// interpolate kd with p = ap1 + bp2 + yp3;
	Color kd(
		a * triangle.v[0].color_diffuse[0] + b * triangle.v[1].color_diffuse[0] + y * triangle.v[2].color_diffuse[0],
		a * triangle.v[0].color_diffuse[1] + b * triangle.v[1].color_diffuse[1] + y * triangle.v[2].color_diffuse[1],
		a * triangle.v[0].color_diffuse[2] + b * triangle.v[1].color_diffuse[2] + y * triangle.v[2].color_diffuse[2]
	
	);
	Color ks(
		a * triangle.v[0].color_specular[0] + b * triangle.v[1].color_specular[0] + y * triangle.v[2].color_specular[0],
		a * triangle.v[0].color_specular[1] + b * triangle.v[1].color_specular[1] + y * triangle.v[2].color_specular[1],
		a * triangle.v[0].color_specular[2] + b * triangle.v[1].color_specular[2] + y * triangle.v[2].color_specular[2]

	);
	double shininess = a * triangle.v[0].shininess + b * triangle.v[1].shininess + y * triangle.v[2].shininess;
	Color color;
	color.r = light.color[0] * (kd.r * LdotN) + ks.r * pow(RdotE, shininess);
	color.g = light.color[1] * (kd.g * LdotN) + ks.g * pow(RdotE, shininess);
	color.b = light.color[2] * (kd.b * LdotN) + ks.b * pow(RdotE, shininess);
	//printf("%f, %f, %f", sphereCenter.x, sphereCenter.y, sphereCenter.z);
	//printf("%f, %f, %f", c.r, c.g, c.b);
	return color;
}

// step3, we do phong shading on sphere,
// fisrt we need to check whether current 
// we'll update the nearestDistance, so that if there is a closer distance in triangle, we will update the return color.
Color checkSpheresIntersect(Ray& ray, Vector3& intersectionPoint, int& intersectSphereID, double& nearestDistance) {
	
	Color final_color = Color();
	// first traverse all spheres to see whether the ray hit one of them.
	// if does, we save its and update our intersect ID and distance
	int tempi = -1;
	float currdis = nearestDistance;
	Vector3 tempPoint = Vector3();
	for (int i = 0; i < num_spheres; i++) {
		Vector3 intersectPoint = Vector3(0.0, 0.0, 0.0);
		// intersectPoint will automatically update because it is a reference..
		if (ray.isIntersectionSphere(spheres[i], intersectPoint) && intersectPoint.z > currdis) {
			tempi = i;
			tempPoint = intersectPoint;
			currdis = intersectPoint.z;
		}
	}
	intersectSphereID = tempi;
	nearestDistance = currdis;
	intersectionPoint = tempPoint;
	if (tempi != -1) {
		final_color = Color();

		for (int j = 0; j < num_lights; j++) {
			Vector3 lightPosition(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
			Vector3 shadowRayDirection = (lightPosition - tempPoint).Normalize();
			// generate a shadow ray to detect intersection;
			Ray shadowRay(tempPoint, shadowRayDirection);
			bool isShadowed = false;
			for (int k = 0; k < num_spheres; k++) {
				Vector3 hitPoint = tempPoint;
				// here should be clear about the shadow ray black and phong model balck...
				// for the shadow ray, it will never hit the same object, since it is a refelct ray from sphere surface..
				// what cause the sphere suface(which opposite to the light source) black is because phong model.
				if (shadowRay.isIntersectionSphere(spheres[k], hitPoint) && k != intersectSphereID) {
					// if the block is behind the light source then there won't be any block
					if (Vector3(lightPosition - tempPoint).len() > Vector3(hitPoint - tempPoint).len()) {
						isShadowed = true;
						break;
					}
				}
			}
			for (int k = 0; k < num_triangles; k++) {
				Vector3 hitPoint = tempPoint;
				if (shadowRay.isIntersectionTriangle(triangles[k], hitPoint)) {
					if (Vector3(lightPosition - tempPoint).len() > Vector3(hitPoint - tempPoint).len()) {
						isShadowed = true;
						break;
					}
				}
			}
			if (!isShadowed) {
				final_color += phongShadingSphere(spheres[tempi], lights[j], tempPoint);
			}
		}

	}

	

	return final_color;

}


Color checkTriangleIntersect(Ray& ray, Vector3& intersectionPoint, int& intersectTriangleID, double& nearestDistance) {

	Color final_color = Color();
	Vector3 intersectPoint = Vector3();
	//int currentTriangleID = -1;
	//int currNearestDistance = nearestDistance;
	// first traverse all triangle to see whether the ray hit one of them.
	// if does, we save its and update our intersect ID and distance
	int tempi = -1;
	float currdis = nearestDistance;
	Vector3 tempPoint = Vector3();
	for (int i = 0; i < num_triangles; i++) {
		// intersectPoint will automatically update because it is a reference..
		if (ray.isIntersectionTriangle(triangles[i], intersectPoint) && intersectPoint.z > currdis) {

			tempi = i;
			tempPoint = intersectPoint;
			currdis = intersectPoint.z;
		}
	}
	intersectTriangleID = tempi;


	intersectionPoint = tempPoint;
	nearestDistance = currdis;

	if (tempi != -1) {
		final_color = Color();
		for (int j = 0; j < num_lights; j++) {
			Vector3 lightPosition(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
			Vector3 shadowRayDirection = (lightPosition - tempPoint).Normalize();
			// generate a shadow ray to detect intersection;
			Ray shadowRay(tempPoint, shadowRayDirection);
			bool isShadowed = false;
			for (int k = 0; k < num_spheres; k++) {
				Vector3 hitPoint;
				// here should be clear about the shadow ray black and phong model balck...
				// for the shadow ray, it will never hit the same object, since it is a refelct ray from sphere surface..
				// what cause the sphere suface(which opposite to the light source) black is because phong model.
				if (shadowRay.isIntersectionSphere(spheres[k], hitPoint)) {
					// if the block is behind the light source then there won't be any block
					if (Vector3(lightPosition - tempPoint).len() > Vector3(hitPoint - tempPoint).len()) {
						isShadowed = true;
						break;
					}
				}
			}
			for (int k = 0; k < num_triangles; k++) {
				Vector3 hitPoint;
				if (shadowRay.isIntersectionTriangle(triangles[k], hitPoint) && (k != tempi)) {
					if (Vector3(lightPosition - tempPoint).len() > Vector3(hitPoint - tempPoint).len()) {
						isShadowed = true;
						break;
					}
				}
			}
			if (!isShadowed) {
				final_color += phongShadingTriangle(triangles[tempi], lights[j], tempPoint);
			}
		}
	}
	
	return final_color;

}

Color computeColor(Ray& r) {
	double nearest = -(double)INT_MAX;
	Vector3 intersect;
	int detectedSphereID = -1, detectedTriangleID = -1;

	Color aggregateColor;
	Color temp;
	temp = checkSpheresIntersect(r, intersect, detectedSphereID, nearest);
	if (detectedSphereID != -1) {
		aggregateColor = temp;
	}
	temp = checkTriangleIntersect(r, intersect, detectedTriangleID, nearest);
	if (detectedTriangleID != -1) {
		aggregateColor = temp;
	}

	if (detectedSphereID == -1 && detectedTriangleID == -1) {
		return Color(1.0, 1.0, 1.0);
	}

	return aggregateColor;
}

int total = 4;
double reflectCoeff = 0.05;
Color computeColorWithReflection(Ray& r, int count) {
	if (count > total) return Color(0.0, 0.0,0.0);
	double nearest = -(double)INT_MAX;
	Vector3 intersect;
	int detectedSphereID = -1, detectedTriangleID = -1;

	Color aggregateColor;
	Color temp;
	temp = checkSpheresIntersect(r, intersect, detectedSphereID, nearest);
	if (detectedSphereID != -1) {
		aggregateColor = temp;
	}
	temp = checkTriangleIntersect(r, intersect, detectedTriangleID, nearest);
	if (detectedTriangleID != -1) {
		aggregateColor = temp;
	}
	//printf("%d", count);
	// if detectedTriangleID != -1, means the cloest object is triangle, so we first calculate it.
	// for the refeclt effection, we need to first the ray direction, and then query the color.
	// for a single reflection, we know N, we Know L , so the relect direction R we can know,
	// the Refect ray will be from intersect point 
	// for the reflection, we only keep tracking for spelular light.
	Ray reflectRay;
	Color ks;
	if (detectedSphereID != -1) {
		// find the sphere, treat light source as a ray shoot from intersect point to the next collider.
		// it just like we ask what will the color be if  shoot a ray from here.
		// what we will get the color, and it create our reflection effect...
		// the reflect ray will be the reflection of our current Ray alone with the normal.
		Sphere sphere = spheres[detectedSphereID];
		Vector3 incident = Vector3() - r.direction;
		Vector3 normal = (intersect - Vector3(sphere.position[0], sphere.position[1], sphere.position[2])).Normalize();
		double NdotL = normal.dot(incident);
		clamp(NdotL);
		Vector3 reflectDir = (normal * (2 * NdotL) - incident).Normalize();
		Vector3 origin = intersect + (reflectDir * SMALL_VALUE);     // add a small value alonge reflectDir so it won't block with current sphere.
		reflectRay = Ray(origin, reflectDir);
		ks = Color(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
	}
	else if (detectedTriangleID != -1) {
		Triangle triangle = triangles[detectedTriangleID];
		Vector3 c0(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
		Vector3 c1(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
		Vector3 c2(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
		Vector3 c = intersect;
		// interpolate normal with  barycentric coordinates method
		Vector3 temp0;
		double areaC0C1C2 = temp0.cross(c2 - c0, c1 - c0).len();
		// a = area(c c1 c2) / area(c0 c1 c2)
		double a = temp0.cross(c2 - c, c1 - c).len() / areaC0C1C2;
		double b = temp0.cross(c2 - c, c0 - c).len() / areaC0C1C2;
		//double y = temp0.cross(c1 - c, c0 - c).len() / areaC0C1C2;
		double y = 1 - a - b;

		// interpolate normal by p = ap1 + bp2 + yp3;
		Vector3 normal = Vector3(
			(a * triangle.v[0].normal[0]) + (b * triangle.v[1].normal[0]) + (y * triangle.v[2].normal[0]),
			(a * triangle.v[0].normal[1]) + (b * triangle.v[1].normal[1]) + (y * triangle.v[2].normal[1]),
			(a * triangle.v[0].normal[2]) + (b * triangle.v[1].normal[2]) + (y * triangle.v[2].normal[2])).Normalize();

		// for the reflection, the light source is 
		Vector3 incident = Vector3() - r.direction;
		double NdotL = normal.dot(incident);
		clamp(NdotL);
		Vector3 reflectDir = (normal * (2 * NdotL) - incident).Normalize();
		Vector3 origin = intersect + (reflectDir * SMALL_VALUE);     // add a small value alonge reflectDir so it won't block with current sphere.
		reflectRay = Ray(origin, reflectDir);
		ks = Color(
			a * triangle.v[0].color_specular[0] + b * triangle.v[1].color_specular[0] + y * triangle.v[2].color_specular[0],
			a * triangle.v[0].color_specular[1] + b * triangle.v[1].color_specular[1] + y * triangle.v[2].color_specular[1],
			a * triangle.v[0].color_specular[2] + b * triangle.v[1].color_specular[2] + y * triangle.v[2].color_specular[2]
		);
	}

	double fact = 0.10;
	if (detectedSphereID == -1 && detectedTriangleID == -1) {
		if (count == 0) {
			return Color(1.0, 1.0, 1.0);
		}
		else {
			return Color(1.0, 1.0, 1.0);
		}
	}
	else {
		Color reflectedResult = computeColorWithReflection(reflectRay, count + 1);
		Color final;
		final.r = aggregateColor.r * (1 - fact) + ks.r * reflectedResult.r *  (fact);
		final.g = aggregateColor.g * (1 - fact) + ks.g * reflectedResult.g *  (fact);
		final.b = aggregateColor.b * (1 - fact) + ks.b * reflectedResult.b *  (fact);
		return final;
	}

	return aggregateColor;
}

// for every light, we simplly remove it to its around, and seperate them to 4 parral light
// to try to create soft effect,
Light* createAreaLight(Light& lights) {
	int nums = min(num_lights, 5);
	Light* areaLight = new Light[nums * 4];
	for (int i = 0; i < nums; i++) {
		// current light[nums]
		areaLight[4 * i];
		areaLight[4 * i];
		areaLight[4 * i];
		areaLight[4 * i + 1];
		areaLight[4 * i + 1];
		areaLight[4 * i + 1];
		areaLight[4 * i + 2];
		areaLight[4 * i + 2];
		areaLight[4 * i + 2];
		areaLight[4 * i + 3];
		areaLight[4 * i + 3];
		areaLight[4 * i + 3];

	}
	return areaLight;
}

// first for the sphere model and without any recursive, we use this model.
Color shootRay(int x, int y) {
	Color finalIntensity;
	if (!enableAntialias) {
		if (!enableReflection) {
			Ray* ray = rayList(x, y);
			finalIntensity = computeColor(ray[0]);
		}
		else{
			Ray* ray = rayList(x, y);
			finalIntensity = computeColorWithReflection(ray[0], 0);
		}
	}
	else {
		if (!enableReflection) {
			Ray* ray = rayList(x, y);
			for (int i = 0; i < 4; i++) {
				finalIntensity += computeColor(ray[i]);
			}
			finalIntensity = finalIntensity / 4;
		}
		else {
			Ray* ray = rayList(x, y);
			for (int i = 0; i < 4; i++) {
				finalIntensity += computeColorWithReflection(ray[i], 0);
			}
			finalIntensity = finalIntensity / 4;
		}
		
	}

	
	return finalIntensity;
}

//MODIFY THIS FUNCTION


// in the drawfunction, we do shoot ray pixel wisely
void draw_scene()
{
	for (unsigned int x = 0; x < WIDTH; x++) {
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (unsigned int y = 0; y < HEIGHT; y++) {
			Color finalIntensity = shootRay(x, y);
			
			// add ambient light
			finalIntensity += Color(ambient_light[0], ambient_light[1], ambient_light[2]);
			finalIntensity.clamp();
			//printf("%f, %f, %f", finalIntensity.r, finalIntensity.g, finalIntensity.b);
			plot_pixel(x, y, finalIntensity.r * 255, finalIntensity.g * 255, finalIntensity.b * 255);
		}
		glEnd();
		glFlush();
	}
	printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc != 5) && (argc != 3) && (argc !=  6) && (argc != 2))
  {  
    printf ("Usage: %s <input scenefile> [enable antialising y] [enable reflction y] [enable softshadow y] [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if (argc == 2) {
	  mode = MODE_DISPLAY;
  }
  else if (argc == 6) {
	  mode = MODE_JPEG;
	  filename = argv[5];
	  char  option = (argv[2])[0];
	  if (option == 'y') {
		  enableAntialias = true;
		  printf("enable antialising! \n");
	  }
	  option = (argv[3])[0];
	  if (option == 'y') {
		  enableReflection = true;
		  printf("enable enableReflection! \n");
	  }

	  option = (argv[4])[0];
	  if (option == 'y') {
		  enableSoftShadow = true;
		  printf("enable enableSoftShadow! \n");
	  }
  }
  else {
	  mode = MODE_DISPLAY;
	  char  option = (argv[2])[0];
	  if (option == 'y') {
		  enableAntialias = true;
		  printf("enable antialising! \n");
	  }
	  option = (argv[3])[0];
	  if (option == 'y') {
		  enableReflection = true;
		  printf("enable enableReflection! \n");
	  }

	  option = (argv[4])[0];
	  if (option == 'y') {
		  enableSoftShadow = true;
		  printf("enable enableSoftShadow! \n");
	  }
  }
  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}


