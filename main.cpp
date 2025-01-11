#include "matrix.hpp"
// #include <iostream>
// #include <vector>
// #include <fstream>
#include <bits/stdc++.h>
using namespace std;

string X_MATRIX_FOLDER = "X_matrix/";
long long int no_of_combinations = 0;
long long int unique_sol_sys = 0;

void print_vector(vector<int> v)
{
    for (auto val : v)
        cout << val << " ";
    cout << endl;
}

void generate_Y_matrix_helper(int n, vector<int> &curr, Matrix<int> &Y)
{
    if (n == 1)
    {
        curr.push_back(1);
        Y.push_back(curr);
        curr.pop_back();
        return;
    }

    curr.push_back(-1);
    generate_Y_matrix_helper(n - 1, curr, Y);
    curr.pop_back();

    curr.push_back(1);
    generate_Y_matrix_helper(n - 1, curr, Y);
    curr.pop_back();

    return;
}

void generate_Y_matrix(int n, Matrix<int> &Y)
{
    vector<int> curr;
    generate_Y_matrix_helper(n, curr, Y);

    cout << "Generated Matrix Y successfully!!" << endl;
    return;
}

void generate_B_matrix(int n, Matrix<int> Y, Matrix<int> &B)
{
    string filename = X_MATRIX_FOLDER + "X_" + to_string(n) + ".txt";
    ifstream file(filename);

    if (!file)
    {
        cerr << "ERROR MESSAGE:: Enable to open the file!!";
        return;
    }

    long long int p = 1;
    string line = "";
    while (getline(file, line) && p < 2 ^ n)
    {
        // x = X[p]
        vector<int> x; // x -> 1 x n
        stringstream ss(line);

        int value;
        while (ss >> value)
        {
            x.push_back(value);
        }

        for (long long int q = 0; q < Y.row(); q++)
        {
            // y = Y[q]
            vector<int> y = Y.getRowVector(q); // y -> 1 x n
            vector<int> b;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    b.push_back(x[j] * y[i]);
                }
            }
            B.push_back(b);
        }
        p++;
    }
    file.close();

    cout << "Generated Matrix B successfully!!" << endl;
    return;
}

void print_solution_matrix(Matrix<int> A, Matrix<long double> t)
{
    std::cout << "[" << std::endl;
    for (unsigned i = 0; i < A.row(); ++i)
    {
        std::cout << "\t[ ";
        for (unsigned j = 0; j < A.col(); ++j)
        {
            std::cout << std::setw(10) << A[i][j] << " ";
        }
        cout << "   |";
        cout << setw(10) << t[i][0] << " ";
        std::cout << std::setw(10) << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
}

// Create augmented matrix: [ A | 1 ]
Matrix<long double> matrix_augmentation(Matrix<int> A)
{
    int A_row = A.row();
    int A_col = A.col();
    Matrix<long double> augmented(A_row, A_col + 1);

    for (size_t i = 0; i < A_row; ++i)
    {
        for (size_t j = 0; j < A_col; ++j)
        {
            augmented[i][j] = static_cast<long double>(A[i][j]);
        }
        augmented[i][A_col] = static_cast<long double>(1);
    }
    return augmented;
}

Matrix<long double> gaussian_elimination(Matrix<int> A)
{
    size_t n = A.row();
    Matrix<long double> augmented = matrix_augmentation(A);

    // Forward Elimination with Partial Pivoting
    for (size_t i = 0; i < n; ++i)
    {
        // Find the pivot row
        size_t pivotRow = i;
        for (size_t k = i + 1; k < n; ++k)
        {
            if (std::abs(augmented[k][i]) > std::abs(augmented[pivotRow][i]))
                pivotRow = k;
        }

        // If the pivot is zero, the matrix is singular
        if (augmented[pivotRow][i] == 0)
        {
            throw std::runtime_error("Matrix is singular or nearly singular!");
        }

        // Swap rows if necessary
        if (pivotRow != i)
        {
            for (size_t j = 0; j <= n; ++j)
                std::swap(augmented[i][j], augmented[pivotRow][j]);
        }

        // Eliminate entries below the pivot
        for (size_t k = i + 1; k < n; ++k)
        {
            long double factor = augmented[k][i] / augmented[i][i];
            for (size_t j = i; j <= n; ++j)
            {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Back Substitution
    Matrix<long double> t(n, 1);
    for (int i = n - 1; i >= 0; --i)
    {
        t[i][0] = augmented[i][n];
        for (size_t j = i + 1; j < n; ++j)
        {
            t[i][0] -= augmented[i][j] * t[j][0];
        }
        t[i][0] /= augmented[i][i];
    }

    return t;
}

// void generate_A_matrix_helper(int n, long long int rowIdx, Matrix<int> B, Matrix<int> A)
// {
//     if (A.row() == n * n)
//     {
//         no_of_combinations++;
//         if (A.determinant() != 0)
//         {
//             unique_sol_sys++;
//             cout << unique_sol_sys << endl;
//             // cout << "Combination of Equations - " << unique_sol_sys << " :: " << endl;
//             // A.printToStdOut();
//             // Gauss-Jordan Elimination::
//             // Matrix<long double> t = gaussian_elimination(A);
//             // print_solution_matrix(A, t);
//         }
//         return;
//     }

//     for (long long int idx = rowIdx; idx < B.row(); idx++)
//     {
//         vector<int> curr_row = B.getRowVector(idx);
//         A.push_back(curr_row);
//         generate_A_matrix_helper(n, idx + 1, B, A);
//         A.pop_back();
//     }
// }

void generate_A_matrix(int n, Matrix<int> B)
{
    no_of_combinations = 0;
    unique_sol_sys = 0;

    vector<int> indices(B.row());
    // Initialize indices with 0, 1, ..., B.row()-1
    iota(indices.begin(), indices.end(), 0);

    // Use a mask to generate combinations of size `n * n`
    vector<bool> mask(B.row(), false);
    fill(mask.begin(), mask.begin() + n * n, true); // First `n * n` entries are true

    // Generate all combinations of size `n*n`
    do
    {
        Matrix<int> A(0, n * n);
        for (size_t i = 0; i < mask.size(); ++i)
        {
            if (mask[i])
            {
                A.push_back(B.getRowVector(indices[i]));
            }
        }

        // A.printToStdOut();
        no_of_combinations++;
        // cout<<no_of_combinations<<endl;
        if (A.determinant() != 0)
        {
            unique_sol_sys++;

            cout << "Unique Solution System #" << unique_sol_sys << ":" << endl;
            // A.printToStdOut();
            Matrix<long double> t = gaussian_elimination(A);
            print_solution_matrix(A, t);
        }
        // for(auto x:mask) cout<<x;
        // cout<<endl;
    } while (prev_permutation(mask.begin(), mask.end()));

    // Matrix<int> A(0, n * n);
    // generate_A_matrix_helper(n, 0, B, A);
    cout << "Generated Matrices A successfully!!" << endl;
}

int main()
{
    try
    {
        long long int no_of_variables = 4;
        cout << "No. of variables in matrix T(size of T):: ";
        cin >> no_of_variables;
        int n = sqrt(no_of_variables); // n*n = no_of_variables

        cout << "Generating matrix Y..." << endl;
        Matrix<int> Y(0, n); // Y -> 2^(n-1) x (n)
        generate_Y_matrix(n, Y);
        cout << "Matrix Y::" << endl;
        Y.printToStdOut();

        cout << "Generating matrix B..." << endl;
        Matrix<int> B(0, no_of_variables); // B -> (2^(2n - 1)) x (n ^ 2)
        generate_B_matrix(n, Y, B);
        cout << "Matrix B::" << endl;
        B.printToStdOut();

        cout << "Generating matrices A..." << endl;
        generate_A_matrix(n, B);
        cout << "No. of combinations(when A might have linearly dependent rows):: " << (1 << (2 * n - 1)) << " C " << B.col() << " = " << no_of_combinations << endl;
        cout << "No. of combinations with unique solution:: " << unique_sol_sys << endl;
    }
    catch (const std::exception &e)
    {
        cerr << "ERROR MESSAGE:: An exception occurred: " << e.what() << endl;
    }

    return 0;
}