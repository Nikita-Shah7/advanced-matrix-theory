#include "matrix.hpp"
#include <bits/stdc++.h>
#include <windows.h>
#include <psapi.h>
#include <omp.h>
using namespace std;

string SOLUTION_SET_FOLDER = "solution_set/";
long long int no_of_combinations = 0;
long long int after_check_1 = 0;
long long int after_check_2 = 0;
long long int after_check_3 = 0;
long long int after_check_4 = 0;
long long int after_check_5 = 0;
long long int non_zero_det_A = 0;
long long int unique_sol_sys = 0;

void printMemoryUsage()
{
    PROCESS_MEMORY_COUNTERS memCounter;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter)))
    {
        std::cout << "Memory usage: " << memCounter.WorkingSetSize / 1024 << " KB\n";
    }
}

long long int nCr(long long int n, long long int r)
{
    if (n < r)
        return 0;
    r = min(r, n - r);
    long long int sum = 1;
    for (int i = 1; i <= r; i++)
        sum = sum * (n - r + i) / i;
    return sum;
}

vector<bool> generate_permutation(int len, int ones, int permutation_idx)
{
    vector<bool> mask(len, false);

    for (size_t i = 0; i < mask.size(); i++)
    {
        if (ones == 0)
            break;

        // Calculate combinations if current position is 1
        long long combinations_with_current_bit = nCr(mask.size() - i - 1, ones - 1);

        if (permutation_idx < combinations_with_current_bit)
        {
            mask[i] = true;
            ones--;
        }
        else
        {
            // Skip to the next combination
            permutation_idx -= combinations_with_current_bit;
        }
    }
    return mask;
}

void print_vector(vector<int> v)
{
    for (auto val : v)
        cout << val << " ";
    cout << endl;
}

template <typename T>
void printToStdOutWithCommas(Matrix<T> A)
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

void generate_B_matrix(int n, const Matrix<int> X, const Matrix<int> Y, Matrix<long double> &B)
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

void print_solution_matrix(Matrix<long double> A, Matrix<long double> t)
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
Matrix<long double> matrix_augmentation(Matrix<long double> A)
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

Matrix<long double> gaussian_elimination(Matrix<long double> A)
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
            // throw std::runtime_error("Matrix is singular or nearly singular!");
            cout << "Matrix is singular or nearly singular!" << endl;
            return {};
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

bool check_condition_for_T(int n, Matrix<int> X, Matrix<long double> t)
{
    if (t.size() == 0)
        return false;

    // Generate matrix T from t
    Matrix<long double> T = generate_T_matrix(n, t);

    for (size_t p = 0; p < X.row(); p++)
    {
        Matrix<long double> x(n, 1); // x -> 1 x n

        for (int k = 0; k < n; k++)
            x[k][0] = X[p][k];

        if ((T * x).l1Norm() > 1)
            return false;
    }
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

bool check_linear_dep_of_A_6_helper(int n, const Matrix<long double> B, vector<bool> mask, int col1, int col2, int col3, int col4)
{
    for (long long int i = 0; i < mask.size(); i++)
    {
        if (mask[i] == 0)
            continue;
        if (B[i][col1] - B[i][col2] != B[i][col3] - B[i][col4])
            return false;
    }
    return true;
}

// Consider columnar dependencies
// C(col1) - C(col2) = C(col3) - C(col4)
bool check_linear_dep_of_A_6(int n, const Matrix<long double> B, vector<bool> mask)
{
    if (n == 2)
        return false;

    if (check_linear_dep_of_A_6_helper(n, B, mask, 0, 1, 3, 4))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 0, 1, 6, 7))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 3, 4, 6, 7))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 1, 2, 4, 5))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 1, 2, 7, 8))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 4, 5, 7, 8))
        return true;
    if (check_linear_dep_of_A_6_helper(n, B, mask, 0, 2, 3, 5))
        return true;

    return false;
}

// check if all entries of a column are same(i.e. all are 1 or -1)
bool check_linear_dep_of_A_7(int n, const Matrix<long double> B, vector<bool> mask)
{
    for (int col = 0; col < B.col(); col++)
    {
        bool is_same = true;
        int val = 0;
        size_t row = 0;
        for (; row < B.row(); row++)
        {
            if (mask[row])
            {
                val = B[row][col];
                break;
            }
        }
        for (; row < B.row(); row++)
        {
            if (mask[row] && B[row][col] != val)
            {
                is_same = false;
                break;
            }
        }
        if (is_same)
            return true;
    }
    return false;
}

void solve_for_A(int n, Matrix<int> X, Matrix<long double> A)
{
    if (A.determinant() != 0)
    {
        non_zero_det_A++;

        // printToStdOutWithCommas(A);
        // A.printToStdOut();
        Matrix<long double> t = gaussian_elimination(A);
        if (check_condition_for_T(n, X, t))
        {
            unique_sol_sys++;
#pragma omp critical
            {
                cout << "Unique Solution System #" << unique_sol_sys << ":" << endl;
                print_solution_matrix(n, t); // file <<
                // print_solution_matrix(A, t); // cout <<
            }
        }
        t.clear();
    }
    return;
}

void generate_A_matrix(int n, Matrix<int> X, Matrix<long double> B)
{
    no_of_combinations = 0;
    after_check_1 = 0;
    after_check_2 = 0;
    after_check_3 = 0;
    after_check_4 = 0;
    after_check_5 = 0;
    non_zero_det_A = 0;
    unique_sol_sys = 0;

    // cout << omp_get_max_threads() << endl;

    omp_set_num_threads(std::min(omp_get_max_threads(), 16));
    omp_set_nested(1);            // Allow nested parallelism
    omp_set_max_active_levels(2); // Limits to 2 levels of parallelism

#ifdef _USING_OMP_
#pragma omp parallel for schedule(dynamic) reduction(+ : no_of_combinations, after_check_1, after_check_2, after_check_3, after_check_4, after_check_5, non_zero_det_A, unique_sol_sys)
#endif
    // for (int ones2 = 0; ones2 <= 0; ones2++)
    for (int ones2 = 1; ones2 <= (n * n) / 2; ones2++)
    {
        int ones1 = n * n - ones2;
        long long int total_permutations1 = nCr(B.row() / 2, ones1);
        long long int total_permutations2 = nCr(B.row() / 2, ones2);

#ifdef _USING_OMP_
#pragma omp parallel for schedule(dynamic) reduction(+ : no_of_combinations, after_check_1, after_check_2, after_check_3, after_check_4, after_check_5, non_zero_det_A, unique_sol_sys)
#endif
        for (long long int permutation1 = 0; permutation1 < total_permutations1; permutation1++)
        {
            // cout << omp_get_thread_num() << endl;
            vector<bool> mask1 = generate_permutation(B.row() / 2, ones1, permutation1);

            for (long long int permutation2 = 0; permutation2 < total_permutations2; permutation2++)
            {
                vector<bool> mask2 = generate_permutation(B.row() / 2, ones2, permutation2);

                vector<bool> mask = mask1;
                mask.insert(mask.end(), mask2.begin(), mask2.end());

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
                Matrix<long double> A(0, n * n);
                for (size_t i = 0; i < mask.size(); ++i)
                {
                    if (mask[i])
                    {
                        A.push_back(B.getRowVector(i));
                    }
                }

                solve_for_A(n, X, A);
                A.clear();
            }
        }
        // cout << ones2 << " " << total_permutations1 << " " << total_permutations2 << endl;
    }
    cout << "Generated Matrices A successfully!!" << endl;
}

int main()
{
    auto start = chrono::high_resolution_clock::now();
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
        // X.printToStdOut();

        cout << "Generating matrix Y..." << endl;
        Matrix<int> Y(static_cast<int>(pow(2, n - 1)), n); // Y -> 2^(n-1) x (n)
        generate_Y_matrix(n, X, Y);
        cout << "Matrix Y::" << endl;
        // Y.printToStdOut();

        cout << "Generating matrix B..." << endl;
        Matrix<long double> B(static_cast<int>(pow(2, 2 * n - 1)), no_of_variables); // B -> (2^(2n - 1)) x (n ^ 2)
        generate_B_matrix(n, X, Y, B);
        cout << "Matrix B::" << endl;
        // B.printToStdOut();
        // printToStdOutWithCommas(B);

        cout << "Generating matrices A..." << endl;
        generate_A_matrix(n, X, B);

        cout << "No. of combinations(when A might have linearly dependent rows):: " << (1 << (2 * n - 1)) << " C " << B.col() << " = " << no_of_combinations << endl;
        cout << "No. of combinations(after check 1):: " << after_check_1 << endl;
        cout << "No. of combinations(after check 2):: " << after_check_2 << endl;
        cout << "No. of combinations(after check 3):: " << after_check_3 << endl;
        cout << "No. of combinations(after check 4):: " << after_check_4 << endl;
        cout << "No. of combinations(after check 5):: " << after_check_5 << endl;
        cout << "No. of combinations(when A has linearly independent rows i.e. det(A)!=0 ):: " << non_zero_det_A << endl;
        cout << "No. of combinations with unique solution:: " << unique_sol_sys << endl;
    }
    catch (const std::exception &e)
    {
        cerr << "ERROR MESSAGE(in main()):: An exception occurred: " << e.what() << endl;
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    cout << "Execution Time: " << duration << " ns" << endl;
    printMemoryUsage();

    return 0;
}