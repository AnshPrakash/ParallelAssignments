#include <iostream>
#include <vector>
#include <cmath>
#include "functions.h"

double EuclidianDistance(Point p1,Point p2){
	return(sqrt(pow(p1.points[0]-p2.points[0],2)+pow(p1.points[1]-p2.points[1],2)+pow(p1.points[2]-p2.points[2],2)));
}
Point add(Point p1,Point p2){
	Point p;
	for (int i = 0; i < sizeof(p1.points)/sizeof(p1.points[0]); ++i){
		p.points[i]=p1.points[i]+p2.points[i];
	}
	return(p);
}
Point divide(Point p1,int k){
	Point p;
	for (int i = 0; i < sizeof(p1.points)/sizeof(p1.points[0]); ++i){
		p.points[i]=p1.points[i]/k;
	}
	return(p);
}
float* vec_to_arr(std::vector<Point> vec){
	float* arr;
	arr=(float*)malloc(sizeof(float)*((vec.size())*3));
	for (int i = 0; i <vec.size(); ++i){
		for (int j = 0; j <3; ++j){
			*(arr+i*3+j)=vec[i].points[j];
		}
	}
	return(arr);
}
std::vector<Point> arr_to_vec(int* arr,int N){
	std::vector<Point> vec;
	for (int i = 0; i <N; ++i){
		Point tmp;
		for (int j = 0;j < 3; ++j){
			tmp.points[j]=*(arr+i*3+j);
		}
		vec.push_back(tmp);
	}
	return(vec);
}