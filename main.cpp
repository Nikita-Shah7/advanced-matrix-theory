#include "matrix.hpp"
// #include <iostream>
// #include <vector>
// #include <fstream>
#include <bits/stdc++.h>
using namespace std;

string X_MATRIX_FOLDER = "X_matrix/";
string SOLUTION_SET_FOLDER = "solution_set/";
long long int no_of_combinations = 0;
long long int after_check_1 = 0;
long long int after_check_2 = 0;
long long int after_check_3 = 0;
long long int after_check_4 = 0;
long long int after_check_5 = 0;
long long int after_check_6 = 0;
long long int non_zero_det_A = 0;
long long int unique_sol_sys = 0;

void print_vector(vector<int> v)
{
    for (auto val : v)
        cout << val << " ";
    cout << endl;
}

void printToStdOutWithCommas(Matrix<int> A)
{
    cout << "[" << endl;
    for (unsigned i = 0; i < A.row(); ++i)
    {
        cout << "\t[ "; // Start of a new row
        for (unsigned j = 0; j < A.col() - 1; ++j)
        {
            cout << setw(1) << A[i][j] << ", ";
        }
        cout << setw(1) << A[i][A.col() - 1] << " ";
        cout << setw(1) << "]," << endl;
    }
    cout << "]" << endl;
}

void generate_X_matrix(int n, Matrix<int> &X)
{
    for (size_t p = 0; p < X.row(); p++)
    {
        X[p][0] = -1;
        for (size_t k = 1; k < n; k++)
        {
            int power = (static_cast<int>(floor(p / pow(2, n - 1 - k)))) % 2;
            X[p][k] = pow(-1, power + 1);
        }
    }
    cout << "Generated HALF Matrix X successfully!!" << endl;
    return;
}

void generate_Y_matrix(int n, const Matrix<int> X, Matrix<int> &Y)
{
    long long int Y_rows = pow(2, n - 1);
    for (size_t q = 0; q < Y_rows; q++)
    {
        for (int l = 0; l < n - 1; l++)
        {
            // int power = (static_cast<int>(floor(q / pow(2, n - 2 - l)))) % 2;
            // Y[q][l] = pow(-1, power + 1);
            Y[q][l] = X[q][l + 1];
        }
        Y[q][n - 1] = 1;
    }

    cout << "Generated Matrix Y successfully!!" << endl;
    return;
}

void generate_B_matrix(int n, const Matrix<int> X, const Matrix<int> Y, Matrix<int> &B)
{
    for (size_t i = 0; i < B.row() / 2; i++)
    {
        int p = static_cast<int>(floor(i / pow(2, n - 1)));
        int q = i % (static_cast<int>(pow(2, n - 1)));
        for (size_t j = 0; j < n * n; j++)
        {
            int k = j % n;
            int l = static_cast<int>(floor(j / n));
            B[i][j] = X[p][k] * Y[q][l];
        }
    }

    for (size_t i = 0; i < B.row() / 2; i++)
    {
        for (size_t j = 0; j < n * n; j++)
        {
            B[i + B.row() / 2][j] = -1 * B[i][j];
        }
    }

    cout << "Generated Matrix B successfully!!" << endl;
    return;
}

void print_solution_matrix(int n, Matrix<long double> t)
{
    string filename = SOLUTION_SET_FOLDER + "sol_set_" + to_string(n) + ".txt";
    // Open the file in append mode
    ofstream file(filename, ios::app);

    if (!file)
    {
        cerr << "ERROR MESSAGE:: Enable to open the file!!";
        return;
    }

    file << "[ ";
    for (int i = 0; i < t.row(); i++)
        file << setw(5) << t[i][0] << " ";
    file << setw(5) << "]" << endl;

    file.close();
}

void print_solution_matrix(Matrix<int> A, Matrix<long double> t)
{
    std::cout << "[" << std::endl;
    for (unsigned i = 0; i < A.row(); ++i)
    {
        std::cout << "\t[ ";
        for (unsigned j = 0; j < A.col(); ++j)
        {
            std::cout << std::setw(5) << A[i][j] << " ";
        }
        cout << "   |";
        cout << setw(5) << t[i][0] << " ";
        std::cout << std::setw(5) << "]" << std::endl;
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

// Generate matrix T from t
Matrix<long double> generate_T_matrix(int n, Matrix<long double> t)
{
    Matrix<long double> T(0, n);

    for (size_t i = 0; i < t.row();)
    {
        vector<long double> v;
        for (size_t j = 0; j < n; j++, i++)
        {
            v.push_back(t[i][0]);
        }
        T.push_back(v);
    }
    return T;
}

bool check_condition_for_T(int n, Matrix<long double> t)
{
    // Generate matrix T from t
    Matrix<long double> T = generate_T_matrix(n, t);

    string filename = X_MATRIX_FOLDER + "X_" + to_string(n) + ".txt";
    ifstream file(filename);

    if (!file)
    {
        cerr << "ERROR MESSAGE(in check condition):: Enable to open the file!!";
        return 0;
    }

    long long int p = 1;
    string line = "";
    while (getline(file, line) && p < pow(2, n - 1))
    {
        // x = X[p]
        Matrix<long double> x(n, 1); // x -> 1 x n
        stringstream ss(line);

        int value;
        int idx = 0;
        while (ss >> value && idx < n)
        {
            x[idx][0] = value;
            idx++;
        }

        if ((T * x).l1Norm() > 1)
            return false;
        p++;
    }
    file.close();
    return true;
}

// In matrix B, ith row = -1*((i+step)th row))
// i.e. these 2 rows are Linearly dep, so these 2 rows can not be taken in same A matrix
bool check_linear_dep_of_A_1(int n, vector<bool> mask)
{
    long long int step = pow(2, 2 * n - 2); // = B.row()/2 = mask.size()/2
    for (long long int i = 0; i < step; i++)
    {
        if (mask[i] && mask[i + step])
            return true;
    }
    return false;
}

bool check_linear_dep_of_A_helper(int idx, vector<bool> mask)
{
    // long long int next_half_step = pow(2, 2 * n - 2); // = B.row()/2 = mask.size()/2
    if (idx >= (mask.size() / 2))
    {
        // cout << idx - (mask.size() / 2) << " ";
        if ((mask[idx] == 1) || (mask[idx - (mask.size() / 2)] == 1))
            return true;
    }
    else
    {
        // cout << idx << " ";
        if ((mask[idx] == 1) || (mask[idx + (mask.size() / 2)] == 1))
            return true;
    }
    return false;
}

// except for n=2 ::
// Ri - R(i+1) = R(i+2) - R(i+3) ; for every i%4=0
// and Ri = -R(i+step) ; for every i
// rows Ri , R(i+1) , R(i+2) , R(i+3) cannot be in same A for every i%4=0
// _______  ________  ________   ________
//    Ri     R(i+1)    R(i+2)     R(i+3)
//    2   *    2    *    2    *    2     =   16 possible cases
// considering n=3 as e.g.
// In matrix B, R0 - R1 = R2 - R3
// rows 0,1,2,3 cannot be in same A,
// rows 0,1,2,19 cannot be in same A,
// rows 0,1,18,19 cannot be in same A,
// rows 16,17,18,3 cannot be in same A,
// etc.
bool check_linear_dep_of_A_2(int n, vector<bool> mask)
{
    if (n == 2)
        return false;

    // long long int next_half_step = pow(2, 2 * n - 2); // = B.row()/2 = mask.size()/2

    for (long long int i = 0; i < mask.size() / 2; i += 4)
    {
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask))
            return true;
    }

    return false;
}

bool check_linear_dep_of_A_3(int n, vector<bool> mask)
{
    if (n == 2)
        return false;

    // long long int next_half_step = pow(2, 2 * n - 2); // = B.row()/2 = mask.size()/2
    long long int step = pow(2, n - 1);

    for (long long int i = 0; i < ((mask.size() / 2) - (3 * step)); i++)
    {
        // p%2=0 should be satisfied where p is for pth row of X
        long long int p = static_cast<int>(floor(((i) / pow(2, n - 1))));
        if (p % 1)
            continue;
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
    }

    return false;
}

// consider linear dependencies with 6 rows (check-4 = check-2 + check-3)
bool check_linear_dep_of_A_4(int n, vector<bool> mask)
{
    if (n == 2)
        return false;

    long long int step = pow(2, n - 1);

    for (long long int i = 0; i < ((mask.size() / 2)); i += 4)
    {
        // p%2=0 should be satisfied for check-3
        long long int p = static_cast<int>(floor(((i) / pow(2, n - 1))));
        if (p % 1)
            continue;
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask))
            return true;
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask))
            return true;
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask))
            return true;
    }

    return false;
}

// consider linear dependencies with 8 rows (check-5 derived from check-4)
bool check_linear_dep_of_A_5(int n, vector<bool> mask)
{
    if (n == 2)
        return false;

    long long int step = pow(2, n - 1);

    for (long long int i = 0; i < ((mask.size() / 2)); i += 4)
    {
        // p%2=0 should be satisfied for check-3
        long long int p = static_cast<int>(floor(((i) / pow(2, n - 1))));
        if (p % 1)
            continue;
        // for i and (i+1)
        if (check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        // for i and (i+2)
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        // for i and (i+3)
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        // for i and (i+step)
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step + 2, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step + 3, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        // for i and (i+2*step)
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step + 2, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step + 3, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step, mask))
            return true;
        // for i and (i+3*step)
        if (check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 * step + 3, mask))
            return true;

        // For (i+1) and (i+1+1)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask))
            return true;
        // For (i+1) and (i+1+2)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask))
            return true;
        // For (i+1) and (i+1+step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step + 2, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask))
            return true;
        // For (i+1) and (i+1+2*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step + 2, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step, mask))
            return true;
        // For (i+1) and (i+1+3*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 1 + 3 * step + 2, mask))
            return true;

        // For (i+2) and (i+2+1)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask))
            return true;
        // For (i+2) and (i+2+step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask))
            return true;
        // For (i+2) and (i+2+2*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step, mask))
            return true;
        // For (i+2) and (i+2+3*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 3, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 2 + 3 * step + 1, mask))
            return true;

        // For (i+3) and (i+3+step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step - 3, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask))
            return true;
        // For (i+3) and (i+3+2*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step - 3, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step - 1, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step, mask))
            return true;
        // For (i+3) and (i+3+3*step)
        if (check_linear_dep_of_A_helper(i, mask) &&
            check_linear_dep_of_A_helper(i + 1, mask) &&
            check_linear_dep_of_A_helper(i + 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 1 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 2 * step, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step - 3, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step - 2, mask) &&
            check_linear_dep_of_A_helper(i + 3 + 3 * step - 1, mask))
            return true;
    }

    return false;
}

void solve_for_A(int n, Matrix<int> A)
{
    if (A.determinant() != 0)
    {
        non_zero_det_A++;

        // printToStdOutWithCommas(A);
        // A.printToStdOut();
        Matrix<long double> t = gaussian_elimination(A);
        if (check_condition_for_T(n, t))
        {
            unique_sol_sys++;
            cout << "Unique Solution System #" << unique_sol_sys << ":" << endl;
            // print_solution_matrix(n, t);  // file <<
            print_solution_matrix(A, t); // cout <<
        }
    }
    return;
}

void generate_A_matrix(int n, Matrix<int> B)
{
    no_of_combinations = 0;
    after_check_1 = 0;
    after_check_2 = 0;
    after_check_3 = 0;
    after_check_4 = 0;
    after_check_5 = 0;
    after_check_6 = 0;
    non_zero_det_A = 0;
    unique_sol_sys = 0;

    // Use a mask to generate combinations of size `n * n`
    vector<bool> mask(B.row(), false);
    fill(mask.begin(), mask.begin() + n * n, true); // First `n * n` entries are true

    // Generate all combinations of size `n*n`
    do
    {
        no_of_combinations++;
        // cout << no_of_combinations << endl;

        if (check_linear_dep_of_A_1(n, mask))
            continue;

        after_check_1++;
        // cout << after_check_1 << endl;

        if (check_linear_dep_of_A_2(n, mask))
            continue;

        after_check_2++;
        // cout << after_check_2 << endl;

        if (check_linear_dep_of_A_3(n, mask))
            continue;

        after_check_3++;
        // cout << after_check_3 << endl;

        if (check_linear_dep_of_A_4(n, mask))
            continue;

        after_check_4++;
        // cout << after_check_4 << endl;

        if (check_linear_dep_of_A_5(n, mask))
            continue;

        after_check_5++;
        // cout << after_check_5 << endl;

        // for (auto x : mask)
        //     cout << x;
        // cout << endl;
        Matrix<int> A(0, n * n);
        for (size_t i = 0; i < mask.size(); ++i)
        {
            if (mask[i])
            {
                A.push_back(B.getRowVector(i));
            }
        }

        // A.printToStdOut();

        // solve_for_A(n, A);

    } while (prev_permutation(mask.begin(), mask.end()));

    cout << "Generated Matrices A successfully!!" << endl;
}

int main()
{
    try
    {
        long long int no_of_variables = 4;
        int n = sqrt(no_of_variables); // n*n = no_of_variables
        do
        {
            cout << "No. of variables in matrix T(size of T):: ";
            cin >> no_of_variables;
            n = sqrt(no_of_variables); // n*n = no_of_variables
        } while (n * n != no_of_variables);

        cout << "Generating HALF matrix X..." << endl;
        Matrix<int> X(static_cast<int>(pow(2, n - 1)), n); // X -> 2^(n) x (n)
        generate_X_matrix(n, X);
        cout << "HALF Matrix X::" << endl;
        X.printToStdOut();

        cout << "Generating matrix Y..." << endl;
        Matrix<int> Y(static_cast<int>(pow(2, n - 1)), n); // Y -> 2^(n-1) x (n)
        generate_Y_matrix(n, X, Y);
        cout << "Matrix Y::" << endl;
        Y.printToStdOut();

        cout << "Generating matrix B..." << endl;
        Matrix<int> B(static_cast<int>(pow(2, 2 * n - 1)), no_of_variables); // B -> (2^(2n - 1)) x (n ^ 2)
        generate_B_matrix(n, X, Y, B);
        cout << "Matrix B::" << endl;
        B.printToStdOut();
        // printToStdOutWithCommas(B);

        cout << "Generating matrices A..." << endl;
        generate_A_matrix(n, B);

        cout << "No. of combinations(when A might have linearly dependent rows):: " << (1 << (2 * n - 1)) << " C " << B.col() << " = " << no_of_combinations << endl;
        cout << "No. of combinations(after check 1):: " << after_check_1 << endl;
        cout << "No. of combinations(after check 2):: " << after_check_2 << endl;
        cout << "No. of combinations(after check 3):: " << after_check_3 << endl;
        cout << "No. of combinations(after check 4):: " << after_check_4 << endl;
        cout << "No. of combinations(after check 5):: " << after_check_5 << endl;
        cout << "No. of combinations(after check 6):: " << after_check_6 << endl;
        cout << "No. of combinations(when A has linearly independent rows i.e. det(A)!=0 ):: " << non_zero_det_A << endl;
        cout << "No. of combinations with unique solution:: " << unique_sol_sys << endl;
    }
    catch (const std::exception &e)
    {
        cerr << "ERROR MESSAGE(in main()):: An exception occurred: " << e.what() << endl;
    }

    return 0;
}