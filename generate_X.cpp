// #include <bits/stdc++.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
using namespace std;

string X_MATRIX_FOLDER = "X_matrix/";
int MAX_n = 4;

void generate_X_matrix_helper(int dim, vector<int> &curr, ofstream &file)
{

    if (dim == 0)
    {
        for (int val : curr)
        {
            file << val << " ";
        }
        file << "\n";
        return;
    }

    curr.push_back(-1);
    generate_X_matrix_helper(dim - 1, curr, file);
    curr.pop_back();

    curr.push_back(1);
    generate_X_matrix_helper(dim - 1, curr, file);
    curr.pop_back();

    return;
}

void generate_X_matrix(int dim)
{
    string filename = X_MATRIX_FOLDER + "X_" + to_string(dim) + ".txt";
    ofstream file(filename);
    if (!file)
    {
        cerr << "ERROR MESSAGE:: Enable to open the file!!";
        return;
    }

    vector<int> curr;
    generate_X_matrix_helper(dim, curr, file);

    file.close();

    cout << "Generated X_" << dim << ".txt successfully!!\n";

    return;
}

int main()
{

    for (int i = 1; i <= MAX_n; i++)
    {
        generate_X_matrix(i);
    }

    return 0;
}