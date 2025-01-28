#include <bits/stdc++.h>
using namespace std;

string SOLUTION_SET_FOLDER = "solution_set/";

// Overload `-` operator for `std::vector<float>`
vector<float> operator-(const vector<float> &vec)
{
    vector<float> negated(vec.size());
    transform(vec.begin(), vec.end(), negated.begin(), [](float x) { 
        if(x==0) return 0.0f;
        return -x; 
    });
    return negated;
}

vector<float> processLine(const string &line)
{
    vector<float> result;
    stringstream ss(line);
    string value;

    while (ss >> value)
    {
        if (value == "[" || value == "]")
            continue;
        float val = stof(value);
        if (val == 0 || val == -0)
        {
            result.push_back(0);
        }
        else
        {
            result.push_back(val);
        }
    }

    return result;
}

void print_solution_set(const set<vector<float>> solution_set)
{
    for (const auto &sol : solution_set)
    {
        cout << "[ ";
        for (auto val : sol)
        {
            cout << setw(7) << val << " ";
        }
        cout << setw(7) << "]" << endl;
    }
}

int main()
{
    long long int no_of_variables = 4;
    int n = sqrt(no_of_variables); // n*n = no_of_variables
    do
    {
        cout << "No. of variables in matrix T(size of T):: ";
        cin >> no_of_variables;
        n = sqrt(no_of_variables); // n*n = no_of_variables
    } while (n * n != no_of_variables);

    // string filename = SOLUTION_SET_FOLDER + "sol_set_" + to_string(n) + ".txt";
    string filename = "unique_sol_3.txt";
    ifstream input_file(filename);
    if (!input_file)
    {
        cerr << "Error: Unable to open file!" << endl;
        return 1;
    }

    long long int no_of_lines = 0;
    set<vector<float>> solution_set;
    string line;

    while (getline(input_file, line))
    {
        vector<float> solution = processLine(line);
        if (solution.empty())
            continue;
        no_of_lines++;
        solution_set.insert(solution);
        solution_set.insert(-solution);
    }

    input_file.close();

    // vector<float> sol(no_of_variables, 0);
    // sol[0] = -1;
    // solution_set.insert(sol);
    // solution_set.insert(-sol);

    // contains half solution only
    print_solution_set(solution_set);
    cout << "No. of Lines Read:: " << no_of_lines << endl;
    cout << "No. of Unique Solutions:: " << solution_set.size() << endl;

    return 0;
}
