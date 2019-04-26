#include "lab3_io.h"
#include "Eigen/Eigen/Dense"

#include <utility>

using namespace Eigen;
void read_matrix (const char* input_filename, int* M, int* N, double** D){
	FILE *fin = fopen(input_filename, "r");
	int i;

	fscanf(fin, "%d%d", M, N);
	
	int num_elements = (*M) * (*N);
	*D = (double*) malloc(sizeof(double)*(num_elements));
	
	for (i = 0; i < num_elements; i++){
		fscanf(fin, "%lf", (*D + i));
	}
	fclose(fin);
}

MatrixXf ptr_to_mat(int M, int N, double *ptr)
{
  MatrixXf ret(M,N);
  for (size_t i = 0; i < M; ++i)
    for (size_t j = 0; j < N; ++j)
      ret(i,j) = *(ptr+i*N+j);
  return ret;
}

std::pair<bool, double> compare_under_tolerance(MatrixXf const& m, MatrixXf const& n, double tolerances[], size_t t_sz)
{
  assert(m.rows() == n.rows());
  assert(m.cols() == n.cols());
  double m_norm = m.norm();
  double n_norm = n.norm();

  double m_norm_abs = fabs(m_norm);
  double n_norm_abs = fabs(n_norm);
  double d = fabs(m_norm_abs-n_norm_abs)*100.0f/n_norm_abs;
  for (size_t i = 0; i < t_sz; ++i) {
    if (d <= tolerances[i]) {
      return std::make_pair(true, tolerances[i]);
    }
  }
  return std::make_pair(false, 100.0f);
}

void write_result (int M, 
		int N, 
		double* D, 
		double* U, 
		double* SIGMA, 
		double* V_T,
		int SIGMAm, 
		int SIGMAn, 
		int K, 
		double* D_HAT,
		double computation_time){
	// Will contain output code
	// printf("hello Cuda\n");
	double tolerances[] = { 0.001, 0.01, 0.1, 1.0 };

  MatrixXf Dm = ptr_to_mat(M, N, D);
  MatrixXf Um = ptr_to_mat(N, N, U);
  MatrixXf Vm = ptr_to_mat(M, M, V_T);
  MatrixXf sigma = MatrixXf::Zero(N,M);
  for (size_t i = 0; i < N; ++i)
    sigma(i, i) = *(SIGMA+i);

  MatrixXf Dm_n = Um*sigma*Vm;
  Dm_n.transposeInPlace();

  auto p = compare_under_tolerance(Dm, Dm_n, tolerances, sizeof(tolerances)/sizeof(double));
  bool d_equal = p.first;
  double d_tolerance = p.second;


  printf("Time: %lf\n", computation_time);
#define FLAG(f) ((f) ? 'T' : 'F')
  printf("%c %lf\n", FLAG(d_equal), d_tolerance);
}
