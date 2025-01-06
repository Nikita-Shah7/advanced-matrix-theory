#include "matrix.hpp"
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#define lli long long int

string X_MATRIX_FOLDER = "X_matrix/";
int no_of_combinations = 0;

void print_vector(vector<int> v)
{
    for (auto val : v)
        cout << val << " ";
    cout << endl;
}

lli nCr(lli n, lli r)
{
    if (n < r)
        return 0;
    lli sum = 1;
    for (int i = 1; i <= r; i++)
        sum = sum * (n - r + i) / i;
    return sum;
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

    cout << "Generated Matrix Y successfully!!\n";
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

    int p = 1;
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

        for (int q = 0; q < Y.row(); q++)
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

    cout << "Generated Matrix B successfully!!\n";
    return;
}

void generate_A_matrix_helper(int n, int Idx, Matrix<int> B, Matrix<int> A)
{
    if (A.row() == n * n)
    {
        no_of_combinations++;
        return;
    }

    for (int idx = Idx; idx < B.row(); idx++)
    {
        A.push_back(B.getRowVector(idx));
        generate_A_matrix_helper(n, idx + 1, B, A);
        A.pop_back();
    }
}

void generate_A_matrix(int n, Matrix<int> B)
{
    Matrix<int> A(0, n * n);
    generate_A_matrix_helper(n, 0, B, A);
}

int main()
{
    try
    {
        int no_of_variables = 4;
        int n = sqrt(no_of_variables); // n*n = no_of_variables

        Matrix<int> Y(0, n); // Y -> 2^(n-1) x (n)
        generate_Y_matrix(n, Y);
        // cout << "Matrix Y::" << endl;
        // Y.printToStdOut();

        Matrix<int> B(0, no_of_variables); // B -> (2^(2n - 1)) x (n ^ 2)
        generate_B_matrix(n, Y, B);
        cout << "Matrix B::" << endl;
        B.printToStdOut();

        generate_A_matrix(n, B);
        cout << "No. of combinations:: " << no_of_combinations << endl;
    }
    catch (const std::exception &e)
    {
        cerr << "ERROR MESSAGE:: An exception occurred: " << e.what() << endl;
    }

    return 0;
}
