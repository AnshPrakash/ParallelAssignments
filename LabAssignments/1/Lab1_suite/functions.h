#pragma once

struct Point{double points[3]={0.0,0.0,0.0};};
double EuclidianDistance(Point p1,Point p2);
Point add(Point p1,Point p2);
Point divide(Point p1,int k);
int* vec_to_arr(std::vector<Point> vec);
std::vector<Point> arr_to_vec(int* arr,int N);