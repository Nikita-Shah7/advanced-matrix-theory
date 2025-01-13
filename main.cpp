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
long long int non_zero_det_A = 0;
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
        cerr << "ERROR MESSAGE(generate_B_matrix):: Enable to open the file!!";
        return;
    }

    long long int p = 0;
    string line = "";
    while (getline(file, line) && p < pow(2, n - 1))
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

    ifstream file2(filename);

    if (!file)
    {
        cerr << "ERROR MESSAGE(generate_B_matrix):: Enable to open the file!!";
        return;
    }

    p = 0;
    line = "";
    while (getline(file2, line) && p < pow(2, n - 1))
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
                    b.push_back(-1 * x[j] * y[i]);
                }
            }
            B.push_back(b);
        }
        p++;
    }
    file2.close();

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
        file << setw(7) << t[i][0] << " ";
    file << setw(10) << "]" << endl;

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
bool check_linear_dep_of_A(int n, vector<bool> mask)
{
    long long int step = pow(2, 2 * n - 2);
    for (long long int i = 0; i < step; i++)
    {
        if (mask[i] && mask[i] == mask[i + step])
            return true;
    }
    return false;
}

void solve_for_A(int n, Matrix<int> A)
{
    if (A.determinant() != 0)
    {
        non_zero_det_A++;

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

        // for (auto x : mask)
        //     cout << x;
        // cout << endl;
        if (check_linear_dep_of_A(n, mask))
            continue;


        Matrix<int> A(0, n * n);
        for (size_t i = 0; i < mask.size(); ++i)
        {
            if (mask[i])
            {
                A.push_back(B.getRowVector(i));
            }
        }

        // A.printToStdOut();
        after_check_1++;
        // cout << after_check_1 << endl;

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
        cout << "No. of combinations(after check 1):: " << after_check_1 << endl;
        cout << "No. of combinations(when A has linearly independent rows i.e. det(A)!=0 ):: " << non_zero_det_A << endl;
        cout << "No. of combinations with unique solution:: " << unique_sol_sys << endl;
    }
    catch (const std::exception &e)
    {
        cerr << "ERROR MESSAGE(in main()):: An exception occurred: " << e.what() << endl;
    }

    return 0;
}