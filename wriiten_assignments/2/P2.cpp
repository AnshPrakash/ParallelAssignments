#include <iostream>
#include <vector>
#include <omp.h>
#include <random>
#include <algorithm>
#include <iterator>
#include <functional>

using namespace std;

void seqsort(std::vector<unsigned int> &X){
	unsigned int i,j,count,N=X.size();
    std::vector<unsigned int> tmp(N);
    for(i=0;i<N;i++){
      count=0;
      for(j=0;j<N;j++)
       if(X[j]<X[i]||X[j]==X[i]&&j<i)
        count++;
      tmp[count]=X[i];
    }
    std::copy(tmp.begin(),tmp.end(),X.begin());
}

void parsort(std::vector<unsigned int> &X){
    unsigned int i,j,count,N=X.size();
    std::vector<unsigned int> tmp(N);

    #pragma omp parallel for schedule(static) shared(tmp) private(j,count) num_threads(2) 
    for(i=0;i<N;i++){
      count=0;
      for(j=0;j<N;j++)
       if(X[j]<X[i]||X[j]==X[i]&&j<i)
        count++;
      tmp[count]=X[i];
    }
    std::copy(tmp.begin(),tmp.end(),X.begin());
}

int main(int argc,char const*argv[]){

	random_device rnd_device;
    // Specify the engine and distribution.
    mt19937 mersenne_engine {rnd_device()};  // Generates random integers
    uniform_int_distribution<int> dist {1, 1000};

    auto gen = [&dist, &mersenne_engine](){
                   return dist(mersenne_engine);
               };

	std::vector< unsigned int> vec(10000);
	generate(begin(vec), end(vec), gen);
	std::vector<unsigned int> vec2=vec;
	for (auto i = vec.begin(); i != vec.end(); ++i)
    	std::cout << *i << ' ';
    std::cout<<"\n";
    double start_time = omp_get_wtime();
	seqsort(vec);
	double time1 = omp_get_wtime() - start_time;
	start_time = omp_get_wtime();
	parsort(vec2);
	double time2 = omp_get_wtime() - start_time;

	for (auto i = vec2.begin(); i != vec2.end(); ++i)
    	std::cout << *i << ' ';
    std::cout<<"\n";
    std::cout<<"Time taken by seqsort "<<time1<<endl;
    std::cout<<"Time taken by parsort "<<time2<<endl;
	return(0);
}