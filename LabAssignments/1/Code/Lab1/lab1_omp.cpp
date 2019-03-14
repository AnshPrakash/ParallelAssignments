/*
*
* Ansh Prakash 2016CS10367
*
*/
#include <iostream>
#include <omp.h> 
#include <vector>
#include <random>
#include "functions.h"
#include "lab1_omp.h"


typedef long long ll;
void kmeans_omp(int num_threads,int N,int K,int* data_points_arr,int** data_point_cluster_arr,float** centroids_arr,int* num_iterations){
	// std::vector<Point> data_points=arr_to_vec(data_points_arr,N);
	std::vector<Point> data_points(N);
	#pragma omp parallel for
	for (int i = 0; i <N; ++i){
		Point tmp;
		for (int j = 0;j < 3; ++j){
			tmp.points[j]=*(data_points_arr+i*3+j);
		}
		data_points[i]=(tmp);
	}
	int n=N; //No. of data points
	int k=K;
	int d=3;
	int nthreads=num_threads;
	int np=1;
	omp_set_num_threads(nthreads);
	std::vector<int> membership(data_points.size());
	std::vector<int> members(k);
	std::vector<Point> centroids(k);
	std::vector<Point> Allcentroids(0);
	int *final_custering;
	final_custering = (int*)malloc(sizeof(int)*((N)*4));
	std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(0,10000);
    //Intialisation of centrioid
	
	for (int i = 0; i < k; ++i){
		Point p;
		for(int j=0;j<d;j++){
			double x=dist6(rng)+double(rand())/(double(RAND_MAX)+1);
			p.points[j]=x;
		}
		centroids[i]=p;
	}
	Allcentroids.insert(std::end(Allcentroids), std::begin(centroids), std::end(centroids));

	// Clustering
	int iters_max=100;
	int iter=0;
	
	while(iter<iters_max){
		iter+=1;
		//assigninng data point to a cluster
		// k is small 
		for (int i = 0; i < k; ++i){
			members[i]=0;
		}
	
		#pragma omp parallel for 
		for (ll i = 0; i < n; ++i){
			double d_min=std::numeric_limits<double>::infinity();
			for (ll j = 0; j < k; ++j){
				double dist=EuclidianDistance(data_points[i],centroids[j]);
				membership[i]=(dist<d_min) ? j:membership[i];
				d_min=(dist<d_min)? dist:d_min;
			}
			#pragma omp atomic
			members[membership[i]]+=1;
		}
		//Recalculating the centroids
		std::vector<Point>(k).swap(centroids);


		// for (ll i = 0; i <n ; ++i){
		// 	centroids[membership[i]]=add(centroids[membership[i]],data_points[i]);
		// }
		std::vector<std::vector<Point>> tp_centroids(nthreads);
		for (int i = 0; i < nthreads; ++i){
			std::vector<Point> temp(1);
			tp_centroids.push_back(temp);
		}
		#pragma omp parallel
		{
			np = omp_get_num_threads();
			int id=omp_get_thread_num();
			std::vector<Point> temp_centroids(k);

			#pragma omp for
			for (ll i = 0; i <n ; ++i){
				temp_centroids[membership[i]]=add(temp_centroids[membership[i]],data_points[i]);
			}

			tp_centroids[id]=temp_centroids;
		}

		for (int i = 0; i <np; ++i){
			for (int j = 0; j < k; ++j){
				centroids[j]=add(centroids[j],tp_centroids[i][j]);
			}
		}

		for (int i = 0; i < k;++i){
			if(members[i]!=0)
				centroids[i]=divide(centroids[i],members[i]);
		}
		Allcentroids.insert(std::end(Allcentroids), std::begin(centroids), std::end(centroids));
	}
	*num_iterations=iters_max;
	// centroids_out (const char* centroid_filename,k,*num_iterations,*centroids);
	for (int i = 0;i<data_points.size(); ++i){
		for (int j= 0; j <d; ++j){
			final_custering[i*4+j]=data_points[i].points[j];
		}
		final_custering[i*4+3]=membership[i];
	}
	*data_point_cluster_arr=final_custering;
	*centroids_arr=vec_to_arr(Allcentroids);
	

}
