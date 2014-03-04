/*
 * main.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: nemo
 */

#include "stdio.h"
#include "math.h"

float rad(float);
float grad(float);

int main(){
	float x;
	float y;
	float z;

	float a;
	float b;
	float b1, b2;
	float c;
	float d;
	float e;
	float f;
	float g;

	float beta;
	float alpha;
	float gamma;
	float delta;
	float fi;
	float epsilon;

	x = 5;
	y = 0;
	z = 0;

	a = 0;
	c = 0;
	beta = rad(10);
	d = 3;
	e = 4;

	b1 = -(a*cos(beta)) + sqrt(a*a*sin(beta)*sin(beta) + x*x + y*y);
	b2 = -(a*cos(beta)) - sqrt(a*a*sin(beta)*sin(beta) + x*x + y*y);

	if(b1>0) b = b1;
	if(b2>0) b = b2;


	delta = atan2(x, y);

	gamma = asin((double)((b*sin(beta))/(sqrt(x*x+y*y))));

	alpha = grad(M_PI/2. - gamma - delta);

	f = fabs(fabs(z) - fabs(c));
	g = sqrt(f*f + b*b);

	gamma = acos((double)((d*d + e*e - g*g)/(2.*d*e)));

	epsilon = asin((double)((e*sin(gamma))/(g)));

	fi = atan2(f, g);

	beta = grad(epsilon - fi);

	gamma = grad(gamma);

	printf("b: %f, alpha: %f, beta: %f, gamma: %f\n", b, alpha, beta, gamma);

	return 0;
}


float rad(float a){
	return (2.*M_PI*a/360.);
}

float grad(float a){
	return (360.*a/(2.*M_PI));
}

