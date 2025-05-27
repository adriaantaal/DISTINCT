#ifndef LEASTSQUARES_H
#define LEASTSQUARES_H

#include <cmath>
#include <iostream>
#include <vector>
#include <numeric> 

#include <fstream>
#include <sstream>
#include <string>

// block size constants
#define N 96
#define LL 3 
#define BLOCK_SIZE 4
#define ADDSIZE 4

// distinct constants
#define maxIter 20
const float EPSILON = 1e-12f; // Small threshold to avoid division by zero
const float STOPCRITERION = 1e-3f; 

// paths
const std::string Hfilename = "../datainputs/H.csv";
const std::string gfilename = "../datainputs/g.csv";
const std::string xtruefilename = "../datainputs/xtrue.csv";
const std::string lambdafilename = "../datainputs/lambda.csv";

// Function declarations
void loadCSVToSquareMat(const std::string& filename, float H[N][N]);
void loadCSVToTallMat(const std::string& filename, float H[N][LL]);
bool loadCSVToLVec(const std::string& filename, float myvec[LL]);
void get_column(const float A[N][LL], float col[N], int j);
void write_column(float A[N][LL], const float col[N], int j);
float frobenius_norm(const float H[N][N]);
float error_norm(const float x[N], const float y[N]);
float compute_xHx(const float H[N][N], const float x[N]);

// Function declarations
void load_block(float A_block[N][BLOCK_SIZE], const float A[N][N], int block_idx);
void generate_test_data(float A[N][N], float b[N], float true_x[N]);
void cholesky_solver(const float H[N][N], const float g[N], float x[N]);

//matrix operations
void compute_AtA_Atb(const float A[N][N], const float b[N], float AtA[N][N], float Atb[N]);
void compute_Atb(const float A[N][N], const float b[N], float Atb[N]);
float vector_dot(const float a[N], const float b[N]);
void copyvec(const float a[N], float b[N]);
void copyintvec(const int a[N], int b[N]);

void merge_unique_indices(const int inds1[N], const int inds2[N], int merged[N], int &merged_K);
void print_NintvecuptoK(const char* label, const int x[N], int K);
void print_Nvector(const char* label, const float x[N]);
void print_LLvector(const char* label, const float x[LL]);
void compute_enabled(const float u[N], bool uenabled[N], int inds[N], int &K, int &KMAX);
void compute_Hxsub(const float Hsub[N][N], const float xsub[N], float Hxsub[N], int K);
void compute_Hx(const float H[N][N], const float x[N], float Hx[N]);
void find_enabled(const float u[N], int inds[N], int &KTOT);
float vectormax(const float U[N]);

// least squares solver, full size
void cholesky_decompose(const float H[N][N], float L[N][N]);
void forward_substitution(const float L[N][N], const float g[N], float y[N]);
void backward_substitution(const float L[N][N], const float y[N], float x[N]);

// least squares solver, upto K size
void cholesky_decompose_K(const float H[N][N], float L[N][N], int K);
void forward_substitution_K(const float L[N][N], const float g[N], float y[N], int K);
void backward_substitution_K(const float L[N][N], const float y[N], float x[N], int K);

void scatter_solution(const float x_sub[N], const int indices[N], int K, float x_full[N]);

float distinct(const float H[N][N], const float g[N], float xout[N], const float yf2, const float lambda);

#endif // LEASTSQUARES_H