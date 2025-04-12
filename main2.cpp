#include "matrix.hpp"
#include <bits/stdc++.h>
#include <windows.h>
#include <psapi.h>
using namespace std;

string SOLUTION_SET_FOLDER = "solution_set/";
long long int no_of_combinations = 0;
long long int non_zero_det_A = 0;
long long int unique_sol_sys = 0;
double FLOATING_POINT_ERROR = 0.001;

void printMemoryUsage()
{
    PROCESS_MEMORY_COUNTERS memCounter;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof(memCounter)))
    {
        std::cout << "Memory usage: " << memCounter.WorkingSetSize / 1024 << " KB\n";
    }
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

bool l1NormCondition(int n, Matrix<long double> T, vector<int> x)
{
    Matrix<long double> Tx(n, 1);

    for (size_t i = 0; i < n; ++i)
    {
        double sum = 0;
        for (size_t k = 0; k < n; ++k)
        {
            sum += T[i][k] * x[k];
        }
        Tx[i][0] = sum;
    }

    return (Tx.l1Norm() <= (1 + FLOATING_POINT_ERROR));
}

bool check_condition_for_T(int n, Matrix<int> X, Matrix<long double> t)
{
    if (t.size() == 0)
        return false;

    // Generate matrix T from t
    Matrix<long double> T = generate_T_matrix(n, t);

    for (size_t p = 0; p < X.row(); p++)
    {
        if (!l1NormCondition(n, T, X.getRowVector(p)))
            return false;
    }
    return true;
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
            cout << "Unique Solution System #" << unique_sol_sys << ":" << endl;
            print_solution_matrix(n, t); // file <<
            // print_solution_matrix(A, t); // cout <<
        }
        t.clear();
    }
    return;
}

bool linear_dependence_relations(int n, const Matrix<long double> mat)
{
    Matrix<long double> mat_transpose = mat.transpose();
    Matrix<long double> nullspace = mat_transpose.kernel();
    // nullspace.printToStdOut();
    if (nullspace.isZero())
        return false;
    return true;
}

bool check_linear_dep_of_X(int n, vector<int> ind, Matrix<int> X)
{
    int A_row = X.row();

    Matrix<long double> Adash(A_row, n);

    for (size_t i = 0; i < A_row; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            Adash[i][j] = static_cast<long double>(X[i][j]);
            Adash[i + A_row][j] = static_cast<long double>(-X[i][j]);
        }
    }
    Matrix<long double> tmp(0, X.col());
    for (auto i : ind)
    {
        tmp.push_back(Adash.getRowVector(i));
    }

    if (linear_dependence_relations(n, tmp))
        return true;
    return false;
}

bool check_linear_dep_of_Y(int n, vector<int> ind, Matrix<int> Y)
{
    int A_row = Y.row();

    Matrix<long double> Adash(A_row, n);

    for (size_t i = 0; i < A_row; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            Adash[i][j] = static_cast<long double>(Y[i][j]);
        }
    }

    Matrix<long double> tmp(0, Y.col());
    for (auto i : ind)
    {
        tmp.push_back(Adash.getRowVector(i));
    }

    if (linear_dependence_relations(n, tmp))
        return true;
    return false;
}

bool check_linear_dep_of_X_1(int n, vector<bool> mask)
{
    for (size_t idx = 0; idx < mask.size() / 2; idx++)
    {
        if ((mask[idx] == 1) && (mask[idx + (mask.size() / 2)] == 1))
            return true;
    }

    return false;
}

bool check_linear_dep_of_X_helper(int idx, vector<bool> mask)
{
    if ((mask[idx] == 1) || (mask[idx + (mask.size() / 2)] == 1))
        return true;

    return false;
}

bool check_linear_dep_of_X_2(int n, vector<bool> mask)
{
    if (n < 4)
        return false;

    for (size_t i = 0; i < mask.size() / 2; i += 4)
    {
        if (check_linear_dep_of_X_helper(i, mask) &&
            check_linear_dep_of_X_helper(i + 1, mask) &&
            check_linear_dep_of_X_helper(i + 2, mask) &&
            check_linear_dep_of_X_helper(i + 3, mask))
            return true;
    }
    return false;
}

bool check_linear_dep_of_X_3(int n, vector<bool> mask)
{
    if (n < 3)
        return false;

    long long int step = pow(2, n - 3);

    for (long long int i = 0; i < ((mask.size() / 2) - (3 * step)); i++)
    {
        // p%2=0 should be satisfied where p is for pth row of X
        // long long int p = static_cast<int>(floor(((i) / pow(2, n - 1))));
        // if (p % 1)
        //     continue;
        if (check_linear_dep_of_X_helper(i, mask) &&
            check_linear_dep_of_X_helper(i + 1 * step, mask) &&
            check_linear_dep_of_X_helper(i + 2 * step, mask) &&
            check_linear_dep_of_X_helper(i + 3 * step, mask))
            return true;
    }
    return false;
}

bool check_linear_dep_of_Y_1(int n, vector<bool> mask)
{
    if (n < 4)
        return false;

    for (size_t i = 0; i < mask.size(); i += 4)
    {
        if (mask[i] &&
            mask[i + 1] &&
            mask[i + 2] &&
            mask[i + 3])
            return true;
    }
    return false;
}

bool check_linear_dep_of_Y_2(int n, vector<bool> mask)
{
    if (n < 3)
        return false;

    long long int step = pow(2, n - 3);

    for (long long int i = 0; i < (mask.size() - (3 * step)); i++)
    {
        // p%2=0 should be satisfied where p is for pth row of X
        // long long int p = static_cast<int>(floor(((i) / pow(2, n - 1))));
        // if (p % 1)
        //     continue;
        if (mask[i] &&
            mask[i + 1 * step] &&
            mask[i + 2 * step] &&
            mask[i + 3 * step])
            return true;
    }
    return false;
}

// generate all possible combination sets from X each with n rows of X
vector<vector<int>> generate_indep_for_X(int n, Matrix<int> X)
{
    vector<vector<int>> x_indep_index;

    int x_rows = pow(2, n);
    long int total_combinations_x = nCr(x_rows, n);

    for (long int combination = 0; combination < total_combinations_x; combination++)
    {
        vector<bool> mask_x = generate_permutation(x_rows, n, combination);

        if (check_linear_dep_of_X_1(n, mask_x))
            continue;

        if (check_linear_dep_of_X_2(n, mask_x))
            continue;

        if (check_linear_dep_of_X_3(n, mask_x))
            continue;

        vector<int> tmp;
        for (long int i = 0; i < mask_x.size(); i++)
        {
            if (mask_x[i])
                tmp.push_back(i);
        }
        if (check_linear_dep_of_X(n, tmp, X))
            continue;

        x_indep_index.push_back(tmp);
    }

    return x_indep_index;
}

// generate all possible combination sets from Y each with n rows of Y
vector<vector<int>> generate_indep_for_Y(int n, Matrix<int> Y)
{
    vector<vector<int>> y_indep_index;

    int y_rows = pow(2, n - 1);
    long int total_combinations_y = nCr(y_rows, n);

    for (long int combination = 0; combination < total_combinations_y; combination++)
    {
        vector<bool> mask_y = generate_permutation(y_rows, n, combination);

        if (check_linear_dep_of_Y_1(n, mask_y))
            continue;

        if (check_linear_dep_of_Y_2(n, mask_y))
            continue;

        vector<int> tmp;
        for (long int i = 0; i < mask_y.size(); i++)
        {
            if (mask_y[i])
                tmp.push_back(i);
        }

        if (check_linear_dep_of_Y(n, tmp, Y))
            continue;

        // // Now generate all permutations of the combination
        // do
        // {
        //     y_indep_index.push_back(tmp);
        // } while (next_permutation(tmp.begin(), tmp.end()));

        y_indep_index.push_back(tmp);
    }
    return y_indep_index;
}

void generate_A_matrix(int n, Matrix<int> X, Matrix<int> Y, Matrix<long double> B)
{
    no_of_combinations = 0;
    non_zero_det_A = 0;
    unique_sol_sys = 0;

    vector<vector<int>> x_indep_index = generate_indep_for_X(n, X);
    vector<vector<int>> y_indep_index = generate_indep_for_Y(n, Y);

    cout << x_indep_index.size() << " " << y_indep_index.size() << endl;

    // long long int total_combi_of_x_per_y_combi = pow(x_indep_index.size(), n);
    // cout << total_combi_of_x_per_y_combi << endl;
    long long int total_combi_of_x_per_y_combi = 0;

    vector<int> x_indep_combi(n, 0);
    while (true)
    {
        // perform code related to current 'x_indep_combi'

        // print_vector(x_indep_combi);
        total_combi_of_x_per_y_combi++;

        for (int y_indep_ind = 0; y_indep_ind < y_indep_index.size(); y_indep_ind++)
        {
            no_of_combinations++;
            // cout << no_of_combinations << endl;

            Matrix<long double> A(0, n * n);
            for (int j = 0; j < n; j++)
            {
                int q = y_indep_index[y_indep_ind][j];
                for (int k = 0; k < n; k++)
                {
                    int p = x_indep_index[x_indep_combi[j]][k];
                    long long int i = p * pow(2, n - 1) + q;
                    // cout << p << " " << q << " " << i << endl;
                    A.push_back(B.getRowVector(i));
                }
            }
            solve_for_A(n, X, A);
            A.clear();
        }

        int pos = n - 1;
        while (pos >= 0 && x_indep_combi[pos] == x_indep_index.size() - 1)
        {
            pos--;
        }

        if (pos < 0)
            break;

        x_indep_combi[pos]++;
        for (int i = pos + 1; i < n; i++)
        {
            x_indep_combi[i] = 0;
        }
    }

    cout << "Generated Matrices A successfully!!" << endl;
    cout << "No. of combinations of X_indep per Y_indep:: " << total_combi_of_x_per_y_combi << endl;
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
        // Y.clear();

        cout << "Generating matrices A..." << endl;
        generate_A_matrix(n, X, Y, B);
        B.clear();

        cout << "No. of combinations(when linearly independent rows of X and Y are considered):: " << no_of_combinations << endl;
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