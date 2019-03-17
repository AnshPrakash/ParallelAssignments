#include <malloc.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

using namespace std;

// Decesending Order
vector<size_t> sort_indices(float* v,int size) {

  // initialize original index locations
  vector<size_t> idx(size);
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

float norm_row_vec(float* v,int n ){
	float norm=0;
	for(int i=0;i<n;i++){
		norm+=v[i]*v[i];
	}
	norm=sqrt(norm);
	return(norm);
}

float inner_product_rv(float* v1,float* v2,int n){
	float inp=0;
	for(int i=0;i< n;i++){
		inp+=v1[i]*v2[i];
	}
	return(inp);
}

void transpose(float* A,float* A_T,int m,int n){
	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j)	{
			A_T[m*j+i]=A[n*i+j];
		}
	}
}



void QRfactorisations(float* A,float* Q,float* R,int n){
	// Assume R is initialised to n X n zero matrix
	// A:dim N X N
	// Assume vectors are represented in col major form in Q,V
	// we will do caculations in row major 
	float* V_T=(float*)(malloc(sizeof(float)*n*n));
	transpose(A,V_T,n,n);
	float* Q_T=(float*)(malloc(sizeof(float)*n*n));
	for (int i=0 ; i<n ; i++){
		R[n*i+i]=norm_row_vec(V_T+n*i,n);//R[i][i]
		for(int j=0;j<n;j++){
			Q_T[n*i+j] = V_T[n*i+j]/R[n*i+i];//R[i][i]
		}
		for (int j=i+1;j<n;j++){
			R[n*i+j] = inner_product_rv(&Q_T[n*i],&V_T[n*j],n);
			for(int k=0;k<n;k++){
				V_T[n*j+k] = V_T[n*j+k] -R[n*i+j]*Q_T[n*i+k];
			}
		}
	}
	transpose(Q_T,Q,n,n);
}

void matmul(float* A,float* B,float* res,int a_m,int a_n,int b_m,int b_n){
	for(int i=0;i<a_m;i++){
		for(int j=0;j<b_n;j++){
			for (int k = 0; k < a_n; ++k){
				res[b_n*i+j] = A[a_n*i+k]*B[b_n*k+j];
			}
		}
	}
}


bool converge(){
	return false;
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T){
	float* SIGMA_p=*SIGMA;
	float* V_T_p=*V_T;
	float* D_T =(float*)(malloc(sizeof(float)*N*M));
	transpose(D,D_T,M,N);
	float* A=(float*)(malloc(sizeof(float)*N*N));
	float* E=(float*)(malloc(sizeof(float)*N*N));
	float* Q=(float*)(malloc(sizeof(float)*N*N));
	float* R=(float*)(malloc(sizeof(float)*N*N));
	matmul(D_T,D,A,N,M,M,N);
	while(!converge()){
		QRfactorisations(A,Q,R,N);
		matmul(R,Q,A,N,N,N,N);
		float* res=(float*)(malloc(sizeof(float)*N*N));
		matmul(E,Q,res,N,N,N,N);
		E=res;
	}
	for (int i = 0; i < N; ++i){
		(*SIGMA)[i]=A[N*i+i];
	}
	float* sig_temp=(float*)(malloc(sizeof(float)*N));
	int idx=0;
	std::vector<size_t> sort_idx = sort_indices((*SIGMA),N);
	for(auto i: sort_idx){
		sig_temp[idx]=(*SIGMA)[i];
		idx+=1;
	}
	for (int i = 0; i < N; ++i){
		(*SIGMA)[i]=sqrt(sig_temp[i]);
	}
	idx=0;
	float* E_T=(float*)(malloc(sizeof(float)*N*N));
	transpose(E,E_T,N,N);
	for (auto i :sort_idx){
		for (int j = 0; j < N; ++j){
			(*U)[N*idx+j]=E_T[N*i+j];
		}
		idx+=1;
	}

	float* sig_inv_ut=(float*)(malloc(sizeof(float)*M*N));//Sigma_inv.U_T
	float* U_T=(float*)(malloc(sizeof(float)*N*N));
	transpose((*U),U_T,N,N);
	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			sig_inv_ut[N*i+j]=U_T[N*i+j]/(*SIGMA)[i];
		}
	}
	for (int i = N; i < M; ++i){
		for (int j = 0; j < N; ++j){
			sig_inv_ut[N*i+j]=0;
		}	
	}
	matmul(sig_inv_ut,D_T,(*V_T),M,N,N,M);
	
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K){
    float* total =(float*)malloc(sizeof(float)*N);
    total[0]=SIGMA[0];
    for (int i = 1; i < N; ++i){
    	total[i]=SIGMA[i]+total[i-1];
    }
    *K=0;
    for (int i = 0; i <N; ++i)    {
    	*K=1+*K;
    	if (total[i]/total[N-1]>= retention){
    		break;
    	}
    }
    *D_HAT=(float*)malloc(sizeof(float)*M*(*K));
    matmul(D,U,(*D_HAT),M,N,N,*K);

}
