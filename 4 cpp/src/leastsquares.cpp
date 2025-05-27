#include "leastsquares.h"


float frobenius_norm(const float H[N][N]) {
    // #pragma HLS INLINE
        float sum = 0.0f;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
    // #pragma HLS PIPELINE
                sum += H[i][j] * H[i][j];
            }
        }
        return sqrtf(sum); // Use sqrtf() for float; synthesizable in HLS
    }

// Load a matrix from a CSV file
void loadCSVToSquareMat(const std::string& filename, float H[N][N]) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open file");

    std::string line;
    int row = 0;
    while (std::getline(file, line) && row < N) {
        std::stringstream ss(line);
        std::string value;
        int col = 0;
        while (std::getline(ss, value, ',') && col < N) {
            H[row][col] = std::stof(value);
            col++;
        }
        row++;
    }
}

void loadCSVToTallMat(const std::string& filename, float H[N][LL]) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open file");

    std::string line;
    int row = 0;
    while (std::getline(file, line) && row < N) {
        std::stringstream ss(line);
        std::string value;
        int col = 0;
        while (std::getline(ss, value, ',') && col < LL) {
            H[row][col] = std::stof(value);
            col++;
        }
        row++;
    }
}

bool loadCSVToLVec(const std::string& filename, float myvec[LL]) {
    std::ifstream file(filename);
    std::string line;
    int i = 0;

    while (std::getline(file, line) && i < LL) {
        if (!line.empty()) {
            myvec[i++] = std::stof(line);
        }
    }

    return (i == N);  // return true only if fully filled
}

// Generate simple test data: A[i][j] = i + j, b[i] = i
void generate_test_data(float A[N][N], float b[N], float true_x[N]) {
    // 1. Generate known solution vector
    for (int i = 0; i < N; ++i) {
        true_x[i] = static_cast<float>((i % 7) + 1); // Any nonzero pattern
    }

    // 2. Create A with linearly independent columns
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            // Ensure independent columns by using a non-repeating base
            A[i][j] = std::sin(0.01f * i * (j + 1)) + 0.1f * (std::rand() % 10);
        }
    }

    // 3. Compute b = A * true_x
    for (int i = 0; i < N; ++i) {
        float sum = 0.0f;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j] * true_x[j];
        }
        b[i] = sum;
    }
}

// Mock memory access for A (you’ll replace this with streaming or BRAM loading)
void load_block(float A_block[N][BLOCK_SIZE], const float A[N][N], int block_idx) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < BLOCK_SIZE; ++j)
            A_block[i][j] = A[i][block_idx * BLOCK_SIZE + j];
}

float vector_dot(const float a[N], const float b[N]) {
    // #pragma HLS INLINE off
        float result = 0.0f;
        for (int i = 0; i < N; ++i) {
    // #pragma HLS PIPELINE II=1
            result += a[i] * b[i];
        }    
        return result;
    }

float error_norm(const float x[N], const float y[N]) {
// #pragma HLS INLINE
    float sum = 0.0f;

    for (int i = 0; i < N; ++i) {
// #pragma HLS PIPELINE
        float diff = x[i] - y[i];
        sum += diff * diff;
    }

    return sqrtf(sum);
}

void compute_AtA_Atb(const float A[N][N], const float b[N], float AtA[N][N], float Atb[N]) {
    // Initialize outputs
    for (int i = 0; i < N; ++i) {
        Atb[i] = 0.0f;
        for (int j = 0; j < N; ++j) {
            AtA[i][j] = 0.0f;
        }
    }

    // Accumulate AtA and Atb
    for (int m = 0; m < N; ++m) {
        for (int i = 0; i < N; ++i) {
            Atb[i] += A[N][i] * b[N];
            for (int j = 0; j < N; ++j) {
                AtA[i][j] += A[N][i] * A[N][j];
            }
        }
    }
}

float compute_xHx(const float H[N][N], const float x[N]) {
    float xHx = 0.0f;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            xHx += x[i] * H[i][j] * x[j];
        }
    }
    return xHx;
}

void compute_Hx(const float H[N][N], const float x[N], float Hx[N]) {
    for (int m = 0; m < N; ++m) {
        for (int i = 0; i < N; ++i) {
            Hx[i] += H[m][i] * x[m];
        }
    }
}

void compute_Atb(const float A[N][N], const float b[N], float Atb[N]) {
    // Initialize outputs
    for (int i = 0; i < N; ++i) {
        Atb[i] = 0.0f;
    }

    // Accumulate Atb
    for (int m = 0; m < N; ++m) {
        for (int i = 0; i < N; ++i) {
            Atb[i] += A[N][i] * b[N];
        }
    }
}

void compute_Hxsub(const float Hsub[N][N], const float xsub[N], float Hxsub[N], int K){
     
    for (int i = 0; i < K; ++i) {
// #pragma HLS PIPELINE II=1
        float acc = 0.0f;
        for (int j = 0; j < K; ++j) {
            acc += Hsub[i][j] * xsub[j];
        }
        Hxsub[i] = acc;
    }
}

void cholesky_solver(const float H[N][N], const float g[N], float x[N]) {
    float L[N][N] = {0};

    // Step 1: Cholesky decomposition: H = L * L^T
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            float sum = H[i][j];
            for (int k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];

            if (i == j)
                L[i][j] = std::sqrt(sum);
            else
                L[i][j] = sum / L[j][j];
        }
    }

    // Step 2: Forward solve L * y = b
    float y[N];
    for (int i = 0; i < N; ++i) {
        float sum = g[i];
        for (int j = 0; j < i; ++j)
            sum -= L[i][j] * y[j];
        y[i] = sum / L[i][i];
    }

    // Step 3: Backward solve L^T * x = y
    for (int i = N - 1; i >= 0; --i) {
        float sum = y[i];
        for (int j = i + 1; j < N; ++j)
            sum -= L[j][i] * x[j];
        x[i] = sum / L[i][i];
    }
}

void cholesky_decompose(const float H[N][N], float L[N][N]) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            L[i][j] = 0.0f;

    for (int k = 0; k < N; ++k) {
        float sum = H[k][k];
        for (int s = 0; s < k; ++s)
            sum -= L[k][s] * L[k][s];
        L[k][k] = std::sqrt(sum);

        for (int i = k + 1; i < N; ++i) {
            float s = H[i][k];
            for (int j = 0; j < k; ++j)
                s -= L[i][j] * L[k][j];
            L[i][k] = s / L[k][k];
        }
    }
}

void forward_substitution(const float L[N][N], const float g[N], float y[N]) {
    // Forward substitution: L y = g
    for (int i = 0; i < N; ++i) {
        float sum = g[i];
        for (int j = 0; j < i; ++j)
            sum -= L[i][j] * y[j];
        y[i] = sum / L[i][i];
    }
}

void backward_substitution(const float L[N][N], const float y[N], float x[N]) {
    // Backward substitution: Lᵗ x = y
    for (int i = N - 1; i >= 0; --i) {
        float sum = y[i];
        for (int j = i + 1; j < N; ++j)
            sum -= L[j][i] * x[j];
        x[i] = sum / L[i][i];
    }
}

void cholesky_decompose_K(const float A[N][N], float L[N][N], int K) {
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j <= i; ++j) {
            float sum = A[i][j];
            for (int k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];

            if (i == j) {
                // Clamp to zero if sum is negative (or just too small)
                L[i][j] = (sum > EPSILON) ? sqrtf(sum) : 0.0f;
            } else {
                // Avoid division by zero
                L[i][j] = (L[j][j] > EPSILON) ? (sum / L[j][j]) : 0.0f;
            }
        }
    }
}

void forward_substitution_K(const float L[N][N], const float b[N], float y[N], int K) {
    for (int i = 0; i < K; ++i) {
        float sum = b[i];
        for (int j = 0; j < i; ++j)
            sum -= L[i][j] * y[j];
        // Avoid division by zero
        y[i] = (L[i][i] > EPSILON) ? (sum / L[i][i]) : 0.0f;
    }
}

void backward_substitution_K(const float L[N][N], const float y[N], float x[N], int K) {
    for (int i = K - 1; i >= 0; --i) {
        float sum = y[i];
        for (int j = i + 1; j < K; ++j)
            sum -= L[j][i] * x[j];
        // Avoid division by zero
        x[i] = (L[i][i] > EPSILON) ? (sum / L[i][i]) : 0.0f;
    }
}

void compute_enabled(const float u[N], bool uenabled[N], int inds[N], int &K, int &KMAX) {
    // #pragma HLS INLINE off

    // Step 1: Compute enabled mask
    for (int i = 0; i < N; ++i) {
            // #pragma HLS UNROLL
        if (u[i] > 0) {
            uenabled[i] = true;
            KMAX += 1;
        } else {
            uenabled[i] = false;
        }
    }

    // Step 2: Initialize top L scores and indices
    float top_vals[ADDSIZE];
    int top_inds[ADDSIZE];
    for (int i = 0; i < ADDSIZE; ++i) {
        top_vals[i] = -1e30f;  // Assume U values are > -1e30
        top_inds[i] = -1;
    }

    // Step 3: Find L largest U[i] with their indices
    for (int i = 0; i < N; ++i) {
        // #pragma HLS PIPELINE II=1
        float val = u[i];

        // Find where to insert in top_vals
        int insert_pos = -1;
        for (int j = 0; j < ADDSIZE; ++j) {
            if (val > top_vals[j]) {
                insert_pos = j;
                break;
            }
        }

        // Shift and insert
        if (insert_pos != -1) {
            for (int j = ADDSIZE - 1; j > insert_pos; --j) {
                top_vals[j] = top_vals[j - 1];
                top_inds[j] = top_inds[j - 1];
            }
            top_vals[insert_pos] = val;
            top_inds[insert_pos] = i;
        }
    }

    // Copy result
    for (int i = 0; i < ADDSIZE; ++i) {
        inds[i] = top_inds[i];
        K = K + 1;
    }
}

void find_enabled(const float u[N], int inds[N], int &KTOT) {
    KTOT = 0;
    int K = 0;
    for (int i = 0; i < N; ++i) {
        if (u[i] > 0) {
            inds[K++] = i;
            KTOT++;
        }
    }
}

void extract_dynamic_submatrix(const float H[N][N], const float g[N], const int indices[N], int K, float H_sub[N][N],  float g_sub[N]){
// #pragma HLS INLINE off
    for (int i = 0; i < K; ++i) {
        g_sub[i] = g[indices[i]];
        for (int j = 0; j < K; ++j) {
// #pragma HLS PIPELINE II=1
            H_sub[i][j] = H[indices[i]][indices[j]];
        }
    }
}

void scatter_solution(const float x_sub[N], const int indices[N], int K, float x_full[N]) {
    for (int i = 0; i < K; ++i) {
// #pragma HLS PIPELINE II=1
        x_full[indices[i]] = x_sub[i];
    }
}

void copyvec(const float a[N], float b[N]){
    for (int i = 0; i < N; ++i) {
        b[i] = a[i];
    }
}

void copyintvec(const int a[N], int b[N]){
    for (int i = 0; i < N; ++i) {
        b[i] = a[i];
    }
}

void sumvec(const float a[N], float b[N]){
    // implements b = b + a
    for (int i = 0; i < N; ++i) {
        b[i] += a[i];
    }
}

void subvec(const float a[N], float b[N]){
    // implements b = b - a
    for (int i = 0; i < N; ++i) {
        b[i] -= a[i];
    }
}

void sumscalar(const float a, float b[N]){
    // implements b = b + a
    for (int i = 0; i < N; ++i) {
        b[i] += a;
    }
}

void convert_int_to_float(const int in[N], float out[N]) {
    // #pragma HLS INLINE off
        for (int i = 0; i < N; ++i) {
    // #pragma HLS UNROLL factor=1 // or use PIPELINE depending on your needs
            out[i] = static_cast<float>(in[i]);
        }
    }

void print_NintvecuptoK(const char* label, const int x[N], int K) {
    std::cout << label << ": [";
    for (int i = 0; i < K; ++i) {
        std::cout << x[i];
        if (i < K - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void print_Nvector(const char* label, const float x[N]) {
    std::cout << label << ": [";
    for (int i = 0; i < N; ++i) {
        std::cout << x[i];
        if (i < N - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void print_LLvector(const char* label, const float x[LL]) {
    std::cout << label << ": [";
    for (int i = 0; i < LL; ++i) {
        std::cout << x[i];
        if (i < LL - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void reset_vec(float x[N]) {
    for (int i = 0; i < N; ++i) {
        x[i] = 0.0f;
    }
}

void reset_mat(float x[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            x[i][j] = 0.0f;
        }
    }
}

void get_column(const float A[N][LL], float col[N], int j) {
    for (int i = 0; i < N; ++i) {
        col[i] = A[i][j];
    }
}

void write_column(float A[N][LL], const float col[N], int j) {
    // #pragma HLS INLINE
        for (int i = 0; i < N; ++i) {
    // #pragma HLS UNROLL
            A[i][j] = col[i];
        }
    }

float vectormax(const float U[N]) {
    float max_val = U[0];
    for (int i = 1; i < N; ++i) {
        if (U[i] > max_val)
            max_val = U[i];
    }   
    return max_val;
}

void merge_unique_indices(const int inds1[N], const int inds2[N],
    int merged[N], int &merged_K
) {
    merged_K = 0;

    // Add all from inds1
    for (int i = 0; i < N; ++i) {
        if (inds1[i] > 0) {  // Optional: use -1 as sentinel if needed
            merged[merged_K++] = inds1[i];
        }
    }

    // Add from inds2 only if not already in merged
    for (int i = 0; i < N; ++i) {
        int val = inds2[i];
        if (val <= 0) continue;  // Optional: skip unused slots

        bool exists = false;
        for (int j = 0; j < merged_K; ++j) {
            if (merged[j] == val) {
                exists = true;
                break;
            }
        }

        if (!exists && merged_K < N) {
            merged[merged_K++] = val;
        }
    }
}

float distinct(const float H[N][N], const float g[N], float xout[N], const float yf2, const float lambda) {
    // Initialize variables
    float L[N][N] = {0};
    float y[N] = {0};
    float u[N] = {0};
    float x[N] = {0};
    float Hx[N] = {0};
    float xsub[N] = {0};

    float Hsub[N][N] = {0};
    float Hxsub[N] = {0};
    float gsub[N] = {0};

    int inds[N] = {0};
    int incrowd[N] = {0};
    int oldcrowd[N] = {0};
    bool uenabled[N] = {0};
    int merged_K = 0;

    // vars to track blocks
    int count = 0;
    int active_blocks = 1; //start with 1 block
    int KMAX = 0; //total number of u>0
    int K = 0; //number of incrowd elements

    // vars to update usefullness and cost function
    copyvec(g, u); //copy g into u
    float ugl[N] = {};
    float J = 0;
    float xHx = 0;
    float xg = 0;
    float xu = 0;

    // vars to calc J
    float xlambda = 0;
    float Jprev = yf2;
    float lambdavec[N] = {0}; 
    sumscalar(lambda,lambdavec); //lambdavec is a vector of lambda
    bool dorun = true;
    // print_Nvector("lambdavec vector", lambdavec);


    for (int ii = 0; ii < maxIter; ++ii) {
        if (!dorun){
            std::cout << "Done at iteration " << ii << ", J = " << J << "\n";
            break;
        }
        else
        {    
            //find ADDSIZE usefullness bigger than 0
            //only include upto toAdd elements

            // Step 1: Get indices of u>0 and count them
            compute_enabled(u, uenabled, inds, K, KMAX);
            // active_blocks = std::floor(count/BLOCK_SIZE); //round down to nearest block size
            // LMAX = active_blocks * BLOCK_SIZE;

            //todo: add new crowd to in crowd
            merge_unique_indices(inds, oldcrowd, incrowd, merged_K);
            copyintvec(incrowd, oldcrowd); //copy incrowd to oldcrowd        

            // reset Hsub and gsub to 0 before extracting
            reset_vec(gsub);
            reset_mat(Hsub);
            extract_dynamic_submatrix(H, g, inds, K, Hsub, gsub);
            // print_Nvector("gsub vector", gsub);

            // set L, y, xsub, Hx to 0 before solving
            reset_vec(y);
            reset_vec(Hx);
            reset_vec(xsub);
            reset_mat(L);
            
            //L starts to have NAN values
            cholesky_decompose_K(Hsub, L, K); // 1. Cholesky: H = L * L^T
            forward_substitution_K(L, gsub, y, K); // 2. Solve L * y = g
            backward_substitution_K(L, y, xsub, K); // 3. Solve L^T * x = y
            
            
            //scatter back to fullsize
            scatter_solution(xsub, inds, K, x); // scatter solution to fullsize x
            // scatter_solution(Hxsub, inds, K, Hx); // scatter solution to fullsize x
            
            // calc Hx in the subspace
            reset_vec(Hx);
            compute_Hx(H, x, Hx); // Hx = H * x
            //todo: do i need to explicitly set Hx non-inds to 0?
            
            //calc new usefulness u = g - H*x   
            reset_vec(u);
            copyvec(g, u); //copy g -> u, i.e. u = g
            subvec(Hx,u); // u = u - H*x

            // print_Nvector("updated u vector", u);
            // print_Nvector("updated x solution", x);

            // calc cost function: J = yf2 - x'*(u + g + lambda)
            xlambda = 0;
            xlambda = vector_dot(x, lambdavec);
            xg = vector_dot(x, g);
            xu = vector_dot(x, u);

            J = yf2 - xu - xg - xlambda;
            std::cout << "Iteration " << ii << ", J = " << J << "\n";

            if ((J < STOPCRITERION) || (std::abs(J - Jprev) < STOPCRITERION)) {dorun = false;}
            Jprev = J;
        }
    }

    copyvec(x, xout); //copy g into u

    return J;
}