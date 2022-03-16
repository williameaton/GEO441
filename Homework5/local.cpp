#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
using namespace std; 


// Aliases 
 using matrix_2d = vector<vector<double>>;

// Function declarations (all beneath main):
double get_k_el_ab(int a, int b, double dx);                              // Calculates k_{ab} for single element
matrix_2d slice_matrix(matrix_2d M, int rows, int cols);                  // Slices rows/cols from end of matrix
matrix_2d matmult(matrix_2d A, matrix_2d B);                              // Matrix multiplication for vector type 

// Functions for inverting a matrix (I shamelessly took these from Geeks4Geeks but edited them a bit to be compatable
// with vectors so that I can use VLAs)
matrix_2d getCofactor(matrix_2d M, matrix_2d temp, int p, int q, int n);
double determinant(matrix_2d M, int n);
matrix_2d adjoint(matrix_2d Q, matrix_2d adj, int dim);
matrix_2d inverse(matrix_2d A, matrix_2d inverse, int N);

// Print functions for debugging
void print_array(double *v, int len);               
void print_2darray(double M[2][2], int row, int col);
void print_vector(vector<double> M);
void print_2dvector(matrix_2d M);



int main(){
    // INITIALISING PARAMETERS:
    // Defining the domain parameters 
    string out_fname = "data.txt";                      // Output file name
    double force = 1.0;                                 // Force (homogenous for this HW hence a scalar)
    int num_el  = 4;                                    // Number of elements
    double q0   = 1.0;                                  // Flux        - LHS boundary condition
    double T1   = 2.0;                                  // Temperature - RHS boundary condition

    double xmin = 0.0;                                  // Domain min x
    double xmax = 1.0;                                  // Domain max x
    double dx   = (xmax-xmin)/double(num_el);           // Grid spacing


    // Create the local shape functions from zeta = -1 to +1 
    vector<double> N1   = {1.0, 0.0};      
    vector<double> N2   = {0.0, 1.0};
    double rhs [2] = {0.0, 0.0};                        // RHS value (Na * F) of vector equation


    double k_el [2][2];                                 // 2d matrix for elemental level stiffness matrix 


    matrix_2d F(num_el+1, vector<double>(1));           // Global scale RHS (force)
    matrix_2d K (num_el+1, vector<double>(num_el+1));   // Gloval scale stiffness matrix
    matrix_2d INVERSE (num_el, vector<double>(num_el)); // Inverse of stiffness matrix

    // **********************************************************************************************************

    // Looping through each element: 
    for (int n=0; n<num_el; n++){

        // Construct local k_ab matrix: 
        for (int a=0; a<2; a++){
            for (int b=0; b<2; b++){
                // Calculate ke_
               k_el[a][b] = get_k_el_ab(a,b, dx);
            }

            // Local RHS which is the integral of Na * f multiplied by dx/2
            rhs[a] = (N1[a] + N2[a])*force * dx/2.0;
        }           

        // Assembly: 
        for (int i=0; i<2; i++){
                for (int j=0; j<2; j++){
                    K[n+i][n+j] += k_el[i][j];
                }
            F[n+i][0] += rhs[i];
        }
    }

    // Add boundary conditions: 
    F[0][0] += q0;
    F[num_el-1][0] += -k_el[0][1]*T1; 

    

    K = slice_matrix(K, 1, 1); // Need to slice out last row and column (n+1) to ensure that matrix can be inverted - this is for the n+1 term that is removed 
    F = slice_matrix(F, 1, 0); // Need to slice out last element 


    // Now need to invert 2D K matrix: 
    INVERSE = inverse(K, INVERSE, num_el);


    // Multiply the K_inv by F to get D (though reusing the F array for memory): 
    F = matmult(INVERSE, F)  ;  


    // Write data to output: 
    // Note that because each of the shape functions are 0 apart from at single points that correspond
    // to the locations of points along the x-array, the inverted d vector is the same as the output U:
    std::ofstream f(out_fname);
    for(int i=0; i<num_el; i++) {
    f << F[i][0] << '\n';
    }

    return 0;
}



















// EXTRA FUNCTIONS: 

double get_k_el_ab(int a, int b, double dx){
    /* k_el_ab is the element of the k_el matrix (elemental stiffness matrix)
     * where a and b can be 1 or 2: 
     * $ k_{ab} = int_{-1}^{+1}[ d_z Na d_z Nb dz/dx  dz]
     * N_a = 1/2 * (1 + ((-1)**a)z)
     * therefore the integral is equal to: 
     * 1/4 x 2 x (-1)^(a+b) x 2/dx where the 2 comes from integrating from -1 to +1 and the 2/dx is realted the the jacobian/change of basis
     */

    return ((1.0/dx) * pow(-1.0, double(a)+double(b)));
}

void print_array(double *v, int len){
    cout << endl; 
    for (int i=0; i<len; i++){
        cout << *(v+i) << ", "; 
    }
    cout << "\n"; 
}

void print_2darray(double M[2][2], int row, int col){

    for (int i=0; i<row; i++){
        for (int j=0; j<col; j++){
            cout << M[i][j] << "  "; 
        }
        cout << "\n";
    }
    cout << "\n"; 
}

void print_vector(vector<double> M){
    cout << endl; 

    for (int j = 0; j < M.size(); j++)
            {
                cout << " " << M[j] << ", ";
            }
            cout << endl;
}

void print_2dvector(matrix_2d M){
    for (int i = 0; i < M.size(); i++)
    {
        for (int j = 0; j < M[i].size(); j++)
        {
            cout << " " << M[i][j] << ", ";
        }
        cout << endl;
    }
}


matrix_2d slice_matrix(matrix_2d M, int rows, int cols){

        int dim_rows = M.size() - rows;
        int dim_cols = M[0].size() - cols;

        matrix_2d P (dim_rows, vector<double>(dim_cols));  

        for (int i = 0; i < dim_rows; i++){
            for (int j = 0; j < dim_cols; j++){
                
                P[i][j] = M[i][j];
            }
        }
    return P;
    }



 // Okay so I shamelessly stole this from Geeks4Geeks as its a nice matrix inverter but modified it for using the vector class:  

// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
matrix_2d getCofactor(matrix_2d M, matrix_2d temp, int p, int q, int n){
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = M[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
}
 



/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
double determinant(matrix_2d M, int n){

    double D = 0.0;        // Initialize result
 
    
    //  Base case : if matrix contains single element
    if (n == 1){
        return M[0][0];}
 
    matrix_2d temp (n-1, vector<double>(n-1));  
 

    double sign = 1.0;  // To store sign multiplier
 
     // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        temp = getCofactor(M, temp, 0, f, n);


        D += sign * M[0][f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint of A[N][N] in adj[N][N].
matrix_2d adjoint(matrix_2d Q, matrix_2d adj, int dim){
    if (dim == 1)
    {
        adj[0][0] = 1;
        return adj;
    }
 
    // temp is used to store cofactors of A[][]
    int sign = 1;
    matrix_2d TEMP (dim, vector<double>(dim));  
 
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<dim; j++)
        {
            // Get cofactor of A[i][j]
            TEMP = getCofactor(Q, TEMP, i, j, dim);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = double((sign))*(determinant(TEMP, dim-1));
        }
    }
    return adj;
}
 

matrix_2d inverse(matrix_2d A, matrix_2d inverse, int N){
    // Find determinant of A[][]
    double det = determinant(A, N);
    if (det == 0.0)
    {
        cout << "Singular matrix, can't find its inverse";
    }
 
    // Find adjoint
    matrix_2d adj (N, vector<double>(N));
    adj = adjoint(A, adj, N);
 

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/float(det);
 
    return inverse;
}


matrix_2d matmult(matrix_2d A, matrix_2d B){
    // Check dimensions: 
    int Arow = A.size();
    int Acol = A[0].size();
    int Brow = B.size();
    int Bcol = B[0].size();
    
    if (Acol != Brow){
        cout << "Matrix dimensions are wrong can not be multiplied.";
    }

    // Create output matrix: 
    matrix_2d M (Arow, vector<double>(Bcol) );


    for (int j = 0; j < Arow; j++){
        for (int k = 0; k < Bcol; k++){
            for (int l = 0; l < Acol; l++){
                M[j][k] += A[j][l]*B[l][k]; 
            }
        }
    }
    return M; 
}
