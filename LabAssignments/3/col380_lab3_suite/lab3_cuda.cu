#include <malloc.h>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

using namespace std;

#define BLOCK_SIZE 8

#define TOLERANCE 0.001
#define JACOBI_UPDATE_TOLERANCE 0.001

double *S; //Symmetric matrix (input)
double  *e; //eigenvalues
double *E; //eigenvectors
int  *ind;
bool *changed;
int  state;
int  N;

void print_vectors(double* vec,int size){
    for (int i = 0; i < size; ++i){
        std::cout<<vec[i]<<" ";
    }
    std::cout<<"\n";
}

void pprint_matrix(double* Mat,int M,int N){
    for (int i = 0; i < M; ++i){
        for (int j = 0; j < N; ++j){
            std::cout<<Mat[N*i+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

double* mat_mul(double* A, int Am, int An, 
                 double* B, int Bm, int Bn){
    double *C;
    C = (double*)malloc(__SIZEOF_DOUBLE__*Am*Bn);

    for (int i=0; i<Am; i++){
        for (int j=0; j<Bn; j++){
            C[i*Bn + j] = 0;
            for (int k=0; k<An; k++){
                C[i*Bn+j] += A[i*An+k] * B[k*Bn+j];
            }
        }
    }

    return C;
}

double* mat_transpose(double* A, int Am, int An) {
    double *B;
    B = (double*)malloc(__SIZEOF_DOUBLE__*An*Am);

    for (int i=0; i<Am; i++){
        for (int j=0; j<An; j++){
            B[ j*Am + i] = A[ i*An + j];
        }
    }

    return B;
}

// __global__ 
// void transpose_kernel(double *A, double *B,int Am,int An);

// __global__ 
// void transpose_kernel(double *A, double *B,int Am,int An){
//     int x = blockIdx.x * blockDim.x + threadIdx.x;
//     int y = blockIdx.y * blockDim.y + threadIdx.y;
//     __syncthreads();
//     if (x >= Am || y >= An)
//         return;
//     B[x*Am + y] = A[y*An +x];
// }

double* cuda_transpose(double* A,int Am,int An){
    double* B;
    B = mat_transpose(A,Am,An);
    return(B);

    // double *dA,*dB;
    // double *B;
    // B = (double*)malloc(sizeof(double)*Am*An);
    // cudaMalloc(&dA,sizeof(double)*Am*An);
    // cudaMalloc(&dB,sizeof(double)*An*Am);

    // cudaMemcpy(A,dA,sizeof(double)*Am*An,cudaMemcpyHostToDevice);

    // dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    // dim3 dimGrid((Am + dimBlock.x - 1) / dimBlock.x,(An + dimBlock.y - 1) / dimBlock.y);

    // transpose_kernel<<<dimGrid,dimBlock>>>(dA,dB,Am,An);
    // cudaThreadSynchronize();
    // cudaMemcpy(dB,B,sizeof(double)*An*Am,cudaMemcpyDeviceToHost);
    // cudaFree(dA);
    // cudaFree(dB);
    // return(B);
}




// __global__ 
// void MatMulKernel(double* A, double* B, double* C,int Am,int An,int Bm,int Bn);

// __global__ 
// void MatMulKernel(double* A, double* B, double* C,int Am,int An,int Bm,int Bn) {
//     double Cvalue = 0;
//     int row = blockIdx.y * blockDim.y + threadIdx.y;
//     int col = blockIdx.x * blockDim.x + threadIdx.x;
//     __syncthreads();

//     if(row > Am || col > Bn) return;
//     for (int e = 0; e < An; ++e)
//         Cvalue += A[row * An + e] *B[e * Bn + col];
//     C[row * Bn + col] = Cvalue;
// }

double* cuda_matmul(double* A,double* B,int Am,int An,int Bm,int Bn){
    double* C;
    C = mat_mul(A,Am,An,B,Bm,Bn);
    return(C);
    // double *dA,*dB,*dC;
    // double* C = (double*)malloc(sizeof(double)*Am*Bn);

    // cudaMalloc(&dA,sizeof(double)*Am*An);
    // cudaMalloc(&dB,sizeof(double)*Bm*Bn);
    // cudaMalloc(&dC,sizeof(double)*Am*Bn);
    // cudaMemcpy(A,dA,sizeof(double)*Am*An,cudaMemcpyHostToDevice);
    // cudaMemcpy(B,dB,sizeof(double)*Bm*Bn,cudaMemcpyHostToDevice);


    // dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    // dim3 dimGrid((Bn + dimBlock.x - 1) / dimBlock.x,(Am + dimBlock.y - 1) / dimBlock.y);
    // MatMulKernel<<<dimGrid,dimBlock>>>(dA,dB,dC,Am,An,Bm,Bn);
    // cudaThreadSynchronize();
    // cudaMemcpy(dC,C,sizeof(double)*Am*Bn,cudaMemcpyDeviceToHost);

    // cudaFree(dA);
    // cudaFree(dB);
    // cudaFree(dC);
    // return(C);

}


int maxind(int k) {
    int m = k+1;

    for (int i = k+2; i < N; i++){
        if (fabs(S[k*N+i]) > fabs(S[k*N+m])){
            m = i;
        }
    }

    return m;
}

void update(int k, double t) {
    double ek_prev = e[k];
    e[k] = ek_prev + t;

    if (e[k] < 0) e[k] = 0;

    if (changed[k] && fabs(ek_prev - e[k]) < JACOBI_UPDATE_TOLERANCE) {
        changed[k] = false;
        state = state - 1;
    }
    else if ((! changed[k]) && fabs(ek_prev - e[k]) > JACOBI_UPDATE_TOLERANCE) {
        changed[k] = true;
        state = state + 1;
    }
}

void rotate(int k, int l, int i, int j, double c, double s,bool eigenvectors){
    double mat1_00 , mat1_01;
    double mat1_10 , mat1_11;
    mat1_00 = c;    mat1_01 = -s;
    mat1_10 = s;    mat1_11 = c;

    double mat2_00 , mat2_10;

    if(eigenvectors){
        mat2_00 = E[i*N + k];
        mat2_10 = E[i*N + l];
    }
    else{
        mat2_00 = S[k*N + l];
        mat2_10 = S[i*N + j];   
    }

    double mat3_00;
    double mat3_10;

    mat3_00 = mat1_00*mat2_00 + mat1_01*mat2_10;
    mat3_10 = mat1_10*mat2_00 + mat1_11*mat2_10;

    if (eigenvectors){
        E[i*N + k] = mat3_00;
        E[i*N + l] = mat3_10;
    }
    else{
        S[k*N + l] = mat3_00;
        S[i*N + j] = mat3_10;
    }
    
}




void print_matrix(double* A, int Am, int An) {
    cout << "[";
    for (int i=0; i<Am; i++){
        if (i>0)
            cout<<" ";
        cout<<"[";
        for (int j=0; j<An-1; j++){
            cout << A[i*An+j] << ", ";
        }
        if (i < Am-1)
            cout << A[i*An+An-1] << "]" << endl;
    }
    cout << A[(Am-1)*An+An-1] << "]]" << endl;
}



// void print_vector(double* A, int An) {
//     cout << "[";
//     for(int i=0; i<An-1; i++)
//         cout << A[i] << ",";
//     cout << A[An-1] << "]" << endl;
// }

void init_jacobi() {
    E = (double*)malloc(sizeof(double)*N*N);
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            E[i*N+j] = 0;
        }
        E[i*N+i] = 1;
    }

    state = N;

    e = (double*)malloc(sizeof(double)*N);
    ind = (int*)malloc(__SIZEOF_INT__*N);
    changed = (bool*)malloc(sizeof(bool)*N);

    for (int k=0; k<N; k++){
        ind[k]     = maxind(k);
        e[k]       = S[k*N+k];
        changed[k] = true;
    }
}

void Jacobi(double *input_matrix, int n,double **eigenvalues, double **eigenvectors) {
    N = n;
    S = input_matrix;

    init_jacobi();

    while(state != 0){
        int m = 0;
        for (int k=1; k<N-1; k++){
            if (fabs(S[k*n+ind[k]]) > fabs(S[m*n+ind[m]])){
                m = k;
            }
        }

        int k = m;
        int l = ind[m];
        double p = S[k*n+l];
        double y = (e[l] - e[k]) / 2.0;
        double d = fabs(y) + sqrt(p*p + y*y);
        double r = sqrt(p*p + d*d);
        double c = d / r;
        double s = p / r;
        double t = (p*p) / d;

        if (y < 0.0) { s = -s; t = -t; }

        S[k*n+l] = 0.0;
        update(k, -t);
        update(l, t);

        for (int i=0; i<k; i++)  { rotate(i, k, i, l, c, s, false); }
        for (int i=k+1; i<l; i++){ rotate(k, i, i, l, c, s, false); }
        for (int i=l+1; i<N; i++)  { rotate(k, i, l, i, c, s, false); }

        for (int i=0; i<N; i++){
            rotate(k, l, i, i, c, s, true);
        }

        ind[k] = maxind(k);
        ind[l] = maxind(l);
    }
    
    *eigenvalues = e;
    *eigenvectors = E;
}

///////////////////////////////////////


// Decesending Order
vector<size_t> sort_indices(double* v,int size) {
  // initialize original index locations
  vector<size_t> idx(size);
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

double norm_row_vec(double* v,int n ){
    double norm=0;
    for(int i=0;i<n;i++){
        norm+=v[i]*v[i];
    }
    norm=sqrt(norm);
    return(norm);
}

double inner_product_rv(double* v1,double* v2,int n){
    double inp=0;
    for(int i=0;i< n;i++){
        inp+=v1[i]*v2[i];
    }
    return(inp);
}








// void QRfactorisations(double* A,double* Q,double* R,int n){
//     // Assume R is initialised to n X n zero matrix
//     // A:dim N X N
//     // Assume vectors are represented in col major form in Q,V
//     // we will do caculations in row major 
//     double* V_T=(double*)(malloc(sizeof(double)*n*n));
//     V_T = mat_transpose(A,n,n);
//     // print_matrix(V_T,n,n);
//     double* Q_T=(double*)(malloc(sizeof(double)*n*n));
//     for (int i=0 ; i<n ; i++){  
//         R[n*i+i]=norm_row_vec(V_T+n*i,n);//R[i][i]
//         for(int j=0;j<n;j++){
//             Q_T[n*i+j] = V_T[n*i+j]/R[n*i+i];//R[i][i]
//         }
//         // #pragma omp parallel for
//         for (int j=i+1;j<n;j++){
//             R[n*i+j] = inner_product_rv(Q_T+n*i,V_T+n*j,n);
//             for(int k=0;k<n;k++){
//                 V_T[n*j+k] = V_T[n*j+k] -R[n*i+j]*Q_T[n*i+k];
//             }
//         }
//     }
//     Q = mat_transpose(Q_T,n,n);
//     // std::cout<<"R\n";
//     // print_matrix(R,n,n);
// }











// void check_result(double *D,double* U,double* SIGMA,double* V_T,int M,int N){
//     double* sig_m=(double*)(malloc(sizeof(double)*N*M));
//     for (int i = 0; i < N; ++i) {
//         for (int j = 0; j < M; ++j){
//             sig_m[M*i+j]=0;
//         }
//     }
//     for (int i = 0; i < N; ++i){
//         sig_m[M*i+i]=SIGMA[i];
//     }
//     double* temp;//=(double*)(malloc(sizeof(double)*N*M));
//     temp =cuda_matmul(U,sig_m,N,N,N,M);
//     double* temp2;//=(double*)(malloc(sizeof(double)*N*M));
//     temp2 =cuda_matmul(temp,V_T,N,M,M,M);
//     std::cout<<"res \n";
//     print_matrix(temp2,N,M);
//     double* D_T;//=(double*)malloc(sizeof(double)*M*N);
//     D_T = transpose(D,M,N);
//     std::cout<<"D_T\n";
//     print_matrix(D_T,N,M);
// }

// /*
//  *****************************************************
//      TODO -- You must implement this function
//  *****************************************************
// */
void SVD(int M, int N, double* D, double** U, double** SIGMA,int SIGMAm,int SIGMAn,double** V_T){
    // double* V_T_p=*V_T;
    double* D_T;// =(double*)(malloc(sizeof(double)*N*M));
    D_T = cuda_transpose(D,M,N);
    double* A;//=(double*)(malloc(sizeof(double)*N*N));
    A = cuda_matmul(D_T,D,N,M,M,N);
    double *eigenvalues,*eigenvectors;
    Jacobi(A,N,(double**)&eigenvalues,(double**)&eigenvectors);
    for (int i = 0; i < N; ++i){
        (*SIGMA)[i] = eigenvalues[i];
    }
    
    double* sig_temp=(double*)(malloc(sizeof(double)*N));
    int idx=0;
    std::vector<size_t> sort_idx = sort_indices((*SIGMA),N);
    for(auto i: sort_idx){
        sig_temp[idx]=(*SIGMA)[i];
        idx+=1;
    }

    for (int i = 0; i < N; ++i){
        (*SIGMA)[i]=sqrt(sig_temp[i]);
    }

    double* E_T;//=(double*)(malloc(sizeof(double)*N*N));
    double* U_T=(double*)(malloc(sizeof(double)*N*N));
    
    // pprint_matrix(*SIGMA,1,N);
    
    
    E_T = cuda_transpose(eigenvectors,N,N);
    

    idx=0;
    for (auto i :sort_idx){
        for (int j = 0; j < N; ++j){
            (U_T)[N*idx+j]=E_T[N*i+j];
        }
        idx+=1;
    }

    *U = U_T;

    double* sig_inv_ut=(double*)(malloc(sizeof(double)*M*N));//Sigma_inv.U_T

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

    (*V_T) = cuda_matmul(sig_inv_ut,D_T,M,N,N,M);
    
}

// /*
//  *****************************************************
//      TODO -- You must implement this function
//  *****************************************************
// */
void PCA(int retention, int M, int N, double* D, double* U, double* SIGMA, double** D_HAT, int *K){
    double* total =(double*)malloc(sizeof(double)*N);
    total[0]=SIGMA[0]*SIGMA[0];
    for (int i = 1; i < N; ++i){
        total[i]=SIGMA[i]*SIGMA[i]+total[i-1];
    }
    *K=0;
    for (int i = 0; i <N; ++i){
        *K=1+*K;
        if ((total[i]/total[N-1])*100>= retention){
            break;
        }
    }
    // *D_HAT=(double*)malloc(sizeof(double)*M*(*K));
    double* U_temp=(double*)malloc(sizeof(double)*N*(*K));

    // #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < (*K); ++j){
            U_temp[(*K)*i+j]=U[N*i+j];
        }
    }

    (*D_HAT) = cuda_matmul(D,U_temp,M,N,N,*K);



    // std::cout<<"U\n";
    // print_matrix(U,N,N);
    // std::cout<<"D_HAT\n";
    // print_matrix((*D_HAT),M,*K);
}

// /*
//  *****************************************************
//      TODO -- You must implement this function
//  *****************************************************
// */
void SVD_and_PCA (int M, 
        int N, 
        double* D, 
        double** U,
        double** SIGMA, 
        double** V_T, 
        int* SIGMAm,
        int* SIGMAn,
        double** D_HAT, 
        int *K,
        int retention) {
    // write your code here
    *SIGMA = (double*)malloc(sizeof(double)*N);
    SVD(M, N, D,U, SIGMA,*SIGMAm,*SIGMAn,V_T);
    PCA(retention, M, N, D, *U, *SIGMA, D_HAT,K);
}

