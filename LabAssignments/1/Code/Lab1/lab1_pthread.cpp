/*
*
* Ansh Prakash 2016CS10367
*
*/
#include <iostream>
// #include <omp.h> 
#include <vector>
#include <pthread.h>
#include <random>
#include "functions.h"
#include "lab1_pthread.h"
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

typedef long long ll;

struct assignment_thread_args{
	int tid;
	int num_threads;
	std::vector<Point> *data_points;
	std::vector<int> *membership;
	std::vector<Point> *centroids;
	std::vector<int> *members;
	ll n;
	int k;

};
void* assignment_thread_function(void* args){
	struct assignment_thread_args s=*(struct assignment_thread_args *)args;
	int num_threads=s.num_threads;
	int tid=s.tid;
	std::vector<Point> &data_points=*s.data_points;
	std::vector<int> &membership=*s.membership;
	std::vector<Point> &centroids=*s.centroids;
	std::vector<int> &members=*s.members;
	ll n=s.n;
	int k=s.k;
	ll chunks=n/num_threads;
	ll last_index=(tid+1)*chunks;
	if(tid==num_threads-1){
		last_index=n;
	}
	for (ll i = (tid*chunks); i < last_index; ++i){
		double d_min=std::numeric_limits<double>::infinity();
		for (ll j = 0; j < k; ++j){
			double dist=EuclidianDistance(data_points[i],centroids[j]);
			membership[i]=(dist<d_min) ? j:membership[i];
			d_min=(dist<d_min)? dist:d_min;
		}
		pthread_mutex_lock( &mutex1 );
		members[membership[i]]+=1;
		pthread_mutex_unlock( &mutex1 );
	}
}

struct recal_centroid_args{
	int tid;
	int num_threads;
	std::vector<Point> *data_points;
	std::vector<int> *membership;
	std::vector<Point> *centroids;
	std::vector<int> *members;
	std::vector<std::vector<Point>> *tp_centroids;
	ll n;
	int k;
};

void* recal_centroid_function(void* args){
	struct recal_centroid_args s=*(struct recal_centroid_args *)args;
	int num_threads=s.num_threads;
	int tid=s.tid;
	std::vector<Point> &data_points=*s.data_points;
	std::vector<int> &membership=*s.membership;
	std::vector<Point> &centroids=*s.centroids;
	std::vector<int> &members=*s.members;
	std::vector<std::vector<Point>> &tp_centroids=*s.tp_centroids;
	ll n=s.n;
	int k=s.k;
	ll chunks=n/num_threads;
	ll last_index=(tid+1)*chunks;
	if(tid==num_threads-1){
		last_index=n;
	}


	int id=tid;
	std::vector<Point> temp_centroids(k);

	
	for (ll i = (tid*chunks); i < last_index; ++i){
		temp_centroids[membership[i]]=add(temp_centroids[membership[i]],data_points[i]);
	}

	tp_centroids[id]=temp_centroids;

}

void kmeans_pthread(int num_threads,int N,int K,int* data_points_arr,int** data_point_cluster_arr,float** centroids_arr,int* num_iterations){
	// std::vector<Point> data_points=arr_to_vec(data_points_arr,N);
	std::vector<Point> data_points(N);
	// #pragma omp parallel for
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
	pthread_t thread_id[nthreads];

	


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
		
		struct assignment_thread_args s1[nthreads];
		for(int i=0; i < num_threads; i++)
	  	{
	  		s1[i].tid=i;
			s1[i].num_threads=nthreads;
			s1[i].data_points=&data_points;
			s1[i].membership=&membership;
			s1[i].centroids=&centroids;
			s1[i].members=&members;
			s1[i].n=n;
			s1[i].k=k;

			pthread_create( &thread_id[i], NULL,assignment_thread_function,(void *)&s1[i]);
	 	}
	 	for (int i = 0; i < num_threads; ++i){
	 		pthread_join( thread_id[i], NULL); 
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
		np=num_threads;
		struct  recal_centroid_args s2[nthreads];
		for(int i=0; i < num_threads; i++){
	  		s2[i].tid=i;
			s2[i].num_threads=nthreads;
			s2[i].data_points=&data_points;
			s2[i].membership=&membership;
			s2[i].centroids=&centroids;
			s2[i].members=&members;
			s2[i].tp_centroids=&tp_centroids;
			s2[i].n=n;
			s2[i].k=k;

			pthread_create(&thread_id[i], NULL,recal_centroid_function,(void *)&s2[i]);
	 	}
	 	for (int i = 0; i < num_threads; ++i){
	 		pthread_join( thread_id[i], NULL); 
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
