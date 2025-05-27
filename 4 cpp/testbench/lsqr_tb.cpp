#include "../src/leastsquares.h"


int main() {
    // Example input matrix A (M x N) and vector b (M)
    // float A[N][N];
    // float b[N];

    // float H[N][N];
    // float g[N];

    // float x[N];

    // // generate test data
    // generate_test_data(A, b, x);
    // compute_AtA_Atb(A,b,H,g);

    // // Call the least squares function
    // cholesky_solver(H,g,x);
    // // checked upto here
    // print_vector("Solution x", x);

    // // start block-decomposition testing
    // float L[N][N];
    // float y[N];
    // float x2[N];

    // cholesky_decompose(H, L); // Step 1: Decompose H = L * Lᵗ   
    // forward_substitution(L, g, y); // Step 2: Solve L y = g_in
    // backward_substitution(L, y, x2); // Step 3: Solve Lᵗ x = y
    // print_vector("Solution to blocked x2", x2);

    // //distinct with test data
    // float xout[N];
    // float J;
    // J = distinct(H, g, xout);

    // print_vector("Solution to distinct xout", xout);
    // std::cout << "Final value for J: " << J << "\n";

    
    // distinct with readout data
    // Load test data from CSV file
    float H[N][N] = {0};
    float gmat[N][LL] = {0};
    float xtruemat[N][LL] = {0};
    float mylambda[LL] = {0};
    float xtrue[N] = {0};
    float g[N] = {0};
    float xout[N] = {0};
    float xoutmat[N][LL] = {0};
    float J[LL] = {0};
    float eps[LL] = {0};
    bool temp[N] = {0};
    int xoutinds[N] = {0};
    int xtrueinds[N] = {0};
    int numxout = 0;
    int numxtrue = 0;
    float yf2 = 0;
    int Ktrue = 0;
    int KtrueMAX = 0;
    float lambda = 0;
    
    loadCSVToSquareMat(Hfilename,H);
    loadCSVToTallMat(gfilename,gmat);
    loadCSVToTallMat(xtruefilename,xtruemat);
    loadCSVToLVec(lambdafilename,mylambda);

    for (int l = 0; l < LL; ++l) {
        get_column(gmat, g, l);
        get_column(xtruemat, xtrue, l);

        lambda = mylambda[l];

        yf2 = compute_xHx(H, xtrue);
        J[l] = yf2 - distinct(H, g, xout, yf2, lambda);
        eps[l] = error_norm(xout, xtrue);
        write_column(xoutmat, xout, l);
        find_enabled(xout, xoutinds, numxout);
        find_enabled(xtrue, xtrueinds, numxtrue);
        std::cout << "Xtrue cardinality: " << numxtrue << " elements in xout: " << numxout << "\n";
        print_NintvecuptoK("Found locations : ", xoutinds, numxout);
        print_NintvecuptoK("True locations : ", xtrueinds, numxtrue);
        std::cout << "Done for = " << l << "\n";
    }

    print_LLvector("Error vector eps = ", eps);
    print_LLvector("Cost function vector J = ", J);

    return 0;
}