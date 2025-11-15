// forward_interp_fixed.cpp
// Newton–Gregory forward interpolation tool (fixed runtime glitches)
// Compile: g++ -O2 -std=c++17 -I /path/to/eigen forward_interp_fixed.cpp -o forward_interp

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

// ANSI color codes (fall back to empty if not supported)
const string GREEN  = "\033[32m";
const string YELLOW = "\033[33m";
const string RED    = "\033[31m";
const string CYAN   = "\033[36m";
const string RESET  = "\033[0m";

static inline bool termSupportsColor() {
    const char* t = getenv("TERM");
    if (!t) return false;
    string ts(t);
    return ts != "dumb";
}

// --------------------------- Helper utilities ---------------------------
bool parseDoublesLine(const string &line, vector<double> &out) {
    out.clear();
    stringstream ss(line);
    double v;
    bool any = false;
    while (ss >> v) { out.push_back(v); any = true; }
    return any;
}

// Read .txt file with two lines: x-line, y-line
bool readFromFileTxt(const string &filename, vector<double> &x, vector<double> &y, string &err) {
    err.clear();
    if (filename.size() < 4 || filename.substr(filename.size()-4) != ".txt") {
        err = "File must have .txt extension.";
        return false;
    }
    ifstream fin(filename);
    if (!fin) { err = "Could not open file."; return false; }
    string lx, ly;
    if (!getline(fin, lx) || !getline(fin, ly)) {
        err = "File format error: must contain two lines (x-line and y-line).";
        return false;
    }
    if (!parseDoublesLine(lx, x)) { err = "Invalid numeric data in first line (x-values)."; return false; }
    if (!parseDoublesLine(ly, y)) { err = "Invalid numeric data in second line (y-values)."; return false; }
    if (x.size() != y.size()) { err = "x and y counts do not match in file."; return false; }
    return true;
}

// Build forward-difference table into Eigen::MatrixXd (n x n, triangular used)
MatrixXd buildForwardTable(const vector<double> &y) {
    int n = (int)y.size();
    MatrixXd D = MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) D(i,0) = y[i];
    for (int j = 1; j < n; ++j)
        for (int i = 0; i + j < n; ++i)
            D(i,j) = D(i+1,j-1) - D(i,j-1);
    return D;
}

string makeASCIITableColored(const MatrixXd &D, int n, bool color) {
    stringstream ss;
    if (color && termSupportsColor()) ss << CYAN;
    ss << "\nForward difference table\n";
    ss << "---------------------------------------------------------------\n";
    ss << setw(12) << "Δ^0";
    for (int j = 1; j < n; ++j) ss << setw(12) << ("Δ^" + to_string(j));
    ss << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n - i; ++j)
            ss << setw(12) << fixed << setprecision(5) << D(i,j);
        ss << "\n";
    }
    if (color && termSupportsColor()) ss << RESET;
    return ss.str();
}

string makeASCIITablePlain(const MatrixXd &D, int n) {
    stringstream ss;
    ss << "Forward difference table\n";
    ss << "---------------------------------------------------------------\n";
    ss << setw(12) << "Δ^0";
    for (int j = 1; j < n; ++j) ss << setw(12) << ("Δ^" + to_string(j));
    ss << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n - i; ++j)
            ss << setw(12) << fixed << setprecision(5) << D(i,j);
        ss << "\n";
    }
    return ss.str();
}

// Stability detector
pair<bool,string> detectStability(const MatrixXd &D, int n) {
    int m_check = min(n-1, 20);
    vector<double> maxabs(m_check+1, 0.0);
    for (int j = 0; j <= m_check; ++j) {
        double m = 0.0;
        for (int i = 0; i + j < n; ++i) m = max(m, fabs(D(i,j)));
        maxabs[j] = m;
    }
    double eps = 1e-12;
    bool all_small = true;
    for (double v : maxabs) if (v > eps) { all_small = false; break; }
    if (all_small) return {true, "Differences are all tiny (stable)."};

    double ratio_threshold = 1e2;
    for (int j = 1; j <= m_check; ++j) {
        double prev = maxabs[j-1], cur = maxabs[j];
        if (prev < 1e-300) continue;
        double ratio = cur / prev;
        if (ratio > ratio_threshold) {
            stringstream ss;
            ss << "Unstable: forward differences grow rapidly (order " << j << " ratio ~ "
               << fixed << setprecision(3) << ratio << ").";
            return {false, ss.str()};
        }
    }
    double max_high = *max_element(maxabs.begin(), maxabs.end());
    if (max_high > 1e12) return {false, "Unstable: very large forward differences encountered."};
    return {true, "No signs of explosive forward-difference growth up to checked order."};
}

// Factorial calculation (more stable than tgamma for small integers)
double factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) result *= i;
    return result;
}

double interpolateNewtonForwardFromTable(const MatrixXd &D, const vector<double> &x, double xp) {
    int n = (int)x.size();
    if (n < 2) return D(0,0);
    double h = x[1] - x[0];
    double s = (xp - x[0]) / h;
    double result = D(0,0);
    double term = 1.0;
    for (int k = 1; k < n; ++k) {
        term *= (s - (k - 1));
        double fact = factorial(k);
        result += (term / fact) * D(0,k);
    }
    return result;
}

double interpolateNewtonForwardTruncated(const MatrixXd &D, const vector<double> &x, double xp, int m) {
    int n = (int)x.size();
    if (n < 2) return D(0,0);
    if (m > n) m = n;
    double h = x[1] - x[0];
    double s = (xp - x[0]) / h;
    double result = D(0,0);
    double term = 1.0;
    for (int k = 1; k < m; ++k) {
        term *= (s - (k - 1));
        double fact = factorial(k);
        result += (term / fact) * D(0,k);
    }
    return result;
}

bool writeResultsToFile(const string &fname, const string &value_line, const string &plain_table, string &err) {
    err.clear();
    vector<string> allowed = {".txt", ".csv", ".md"};
    bool okext = false;
    for (auto &e : allowed)
        if (fname.size() > e.size() && fname.substr(fname.size()-e.size()) == e) { okext = true; break; }
    if (!okext) { err = "Output filename must use .txt, .csv or .md extension."; return false; }
    ofstream fout(fname);
    if (!fout) { err = "Cannot open file for writing."; return false; }
    fout << value_line << "\n\n" << plain_table << "\n";
    fout.close();
    return true;
}

// ------------------------------ Main -------------------------------------
int main() {
    // Removed ios::sync_with_stdio(false) and cin.tie(nullptr)
    // to fix prompt display issue when using cout before cin

    bool color = termSupportsColor();
    if (color) cout << CYAN;
    cout << "Newton–Gregory Forward Interpolation Tool\n";
    cout << "-----------------------------------------------------------\n";
    cout << "Input methods:\n";
    cout << " • Direct input (space-separated x-values and y-values)\n";
    cout << " • File input (use when you have many points): file must be .txt with 2 lines\n";
    cout << "   Line 1: space-separated x-values\n";
    cout << "   Line 2: space-separated y-values\n";
    cout << "The program prints the interpolated value to 5 decimal places and can save results to a file.\n";
    if (color) cout << RESET;

    vector<double> x, y;
    string line;

    // Read x-values (direct line). If >10 values, switch to file prompt
    while (true) {
        if (color) cout << GREEN;
        cout << "Enter x-values (space-separated): " << (color ? RESET : "");
        if (!getline(cin, line)) return 0;
        if (!parseDoublesLine(line, x)) {
            if (color) cout << RED;
            cout << "Invalid input. Please re-enter x-values as space-separated numbers.\n" << (color?RESET:"");
            continue;
        }
        break;
    }

    // If user provided >10 via direct typing, switch to file-based input
    if (x.size() > 10) {
        if (color) cout << YELLOW;
        cout << "You entered more than 10 x-values. Please provide a .txt file (two lines: x-line and y-line).\n" << (color?RESET:"");
        string fname;
        while (true) {
            if (color) cout << GREEN;
            cout << "Enter filename (with .txt): " << (color?RESET:"");
            if (!getline(cin, fname)) return 0;
            string ferr;
            if (!readFromFileTxt(fname, x, y, ferr)) {
                if (color) cout << RED;
                cout << "File error: " << ferr << "\n" << (color?RESET:"");
                cout << "Please correct the file and enter filename again.\n";
                continue;
            }
            break;
        }
    } else {
        // Direct mode: get y-values directly (must match count)
        while (true) {
            if (color) cout << GREEN;
            cout << "Enter corresponding y-values (space-separated): " << (color?RESET:"");
            if (!getline(cin, line)) return 0;
            if (!parseDoublesLine(line, y) || y.size() != x.size()) {
                if (color) cout << RED;
                cout << "Invalid input. Ensure the y-values are numeric and the count matches x-values.\n" << (color?RESET:"");
                continue;
            }
            break;
        }
    }

    int n = (int)x.size();

    // If only 1 point — return that value directly (still keep y present)
    if (n == 1) {
        if (y.empty()) {
            if (color) cout << RED;
            cout << "Error: single x provided but no y-value. Please restart and provide both.\n" << (color?RESET:"");
            return 0;
        }
        if (color) cout << YELLOW;
        cout << "Only one data point provided. Interpolated value equals the given y-value.\n" << (color?RESET:"");
    }

    if (y.empty()) {
        if (color) cout << RED;
        cout << "Internal error: y-values missing.\n" << (color?RESET:"");
        return 0;
    }

    // Very large tables warning
    if (n > 5000) {
        if (color) cout << YELLOW;
        cout << "Warning: you provided " << n << " points. Very large tables may be slow and numerically unstable.\n";
        cout << "Do you want to proceed? (yes/no): " << (color?RESET:"");
        string ans;
        if (!getline(cin, ans)) return 0;
        if (ans != "yes" && ans != "Yes" && ans != "y" && ans != "Y") {
            cout << "Operation cancelled by user.\n";
            return 0;
        }
    }

    // Check strictly increasing x
    for (int i = 1; i < n; ++i) {
        if (!(x[i] > x[i-1])) {
            if (color) cout << RED;
            cout << "Error: x-values must be strictly increasing. Please fix your data and try again.\n" << (color?RESET:"");
            return 0;
        }
    }

    if (n >= 2) {
        double h = x[1] - x[0];
        double tol = fabs(h) * 1e-8;
        bool eq = true;
        for (int i = 2; i < n; ++i) {
            if (fabs((x[i] - x[i-1]) - h) > tol) { eq = false; break; }
        }
        if (!eq) {
            if (color) cout << RED;
            cout << "Spacing is not uniform within tolerance. Forward interpolation requires equal spacing.\n" << (color?RESET:"");
            return 0;
        }
    }

    // Read interpolation point xp
    double xp;
    while (true) {
        if (color) cout << GREEN;
        cout << "Enter the point at which you want the interpolated value: " << (color?RESET:"");
        if (!getline(cin, line)) return 0;
        stringstream ss(line);
        if (ss >> xp) break;
        if (color) cout << RED;
        cout << "Invalid numeric input. Please re-enter the point.\n" << (color?RESET:"");
    }

    // Decide full vs truncated table
    int full_table_threshold = 1000;
    bool use_full_table = true;
    if (n > full_table_threshold) {
        if (color) cout << YELLOW;
        cout << "Large dataset (" << n << " points). Building full forward-difference table is O(n^2) memory.\n";
        cout << "This may be slow and use a lot of RAM. Proceed building full table? (yes/no) [no = use truncated]: " << (color?RESET:"");
        string ans;
        if (!getline(cin, ans)) return 0;
        if (ans != "yes" && ans != "Yes" && ans != "y" && ans != "Y") {
            use_full_table = false;
        }
    }

    MatrixXd D;
    int used_n = n;
    if (use_full_table) {
        try {
            D = MatrixXd::Zero(n, n);
            for (int i = 0; i < n; ++i) D(i,0) = y[i];
            for (int j = 1; j < n; ++j)
                for (int i = 0; i + j < n; ++i)
                    D(i,j) = D(i+1,j-1) - D(i,j-1);
            used_n = n;
        } catch (bad_alloc &e) {
            if (color) cout << RED;
            cout << "Memory allocation failed while building full forward-difference table.\n" << (color?RESET:"");
            return 0;
        }
    } else {
        int m = min(n, 200);
        if (color) cout << YELLOW;
        cout << "Using truncated approach with first " << m << " points (approximate result).\n" << (color?RESET:"");
        used_n = m;
        D = MatrixXd::Zero(m, m);
        for (int i = 0; i < m; ++i) D(i,0) = y[i];
        for (int j = 1; j < m; ++j)
            for (int i = 0; i + j < m; ++i)
                D(i,j) = D(i+1,j-1) - D(i,j-1);
    }

    // Stability detection
    auto stab = detectStability(D, used_n);
    if (!stab.first) {
        if (color) cout << YELLOW;
        cout << "Stability warning: " << stab.second << "\n";
        cout << "High-degree Newton polynomial may be numerically unstable. Consider using fewer points or a spline.\n" << (color?RESET:"");
        if (color) cout << GREEN;
        cout << "Do you want to continue with current method? (yes/no): " << (color?RESET:"");
        string ans;
        if (!getline(cin, ans)) return 0;
        if (ans != "yes" && ans != "Yes" && ans != "y" && ans != "Y") {
            cout << "Operation cancelled by user due to instability concerns.\n";
            return 0;
        }
    } else {
        if (color) cout << CYAN;
        cout << "Stability check passed: " << stab.second << "\n" << (color?RESET:"");
    }

    // Interpolation
    double result;
    if (n == 1) result = y[0];
    else if (use_full_table || used_n == n) result = interpolateNewtonForwardFromTable(D, x, xp);
    else result = interpolateNewtonForwardTruncated(D, x, xp, used_n);

    if (color) cout << YELLOW;
    cout << "\nInterpolated value at x = " << fixed << setprecision(5) << xp << " is: " << result << "\n";
    if (color) cout << RESET;

    int print_table_rows = used_n;
    const int TABLE_PRINT_LIMIT = 60;
    if (used_n > TABLE_PRINT_LIMIT) {
        print_table_rows = TABLE_PRINT_LIMIT;
        if (color) cout << YELLOW;
        cout << "(Large table truncated for display; saved output file will contain more.)\n";
        if (color) cout << RESET;
    }
    MatrixXd Dprint = D.topLeftCorner(print_table_rows, print_table_rows);
    cout << makeASCIITableColored(Dprint, print_table_rows, color) << "\n";

    if (color) cout << GREEN;
    cout << "Would you like to save the result and the full difference table to a file? (yes/no): " << (color?RESET:"");
    string saveAns;
    if (!getline(cin, saveAns)) return 0;
    if (saveAns == "yes" || saveAns == "Yes" || saveAns == "y" || saveAns == "Y") {
        string outname;
        while (true) {
            if (color) cout << GREEN;
            cout << "Enter output filename (use .txt, .csv or .md): " << (color?RESET:"");
            if (!getline(cin, outname)) return 0;
            string err;
            string plain_table = makeASCIITablePlain(D, used_n);
            string value_line = "Interpolated value at x = " + to_string(xp) + " is: " + to_string(result);
            if (writeResultsToFile(outname, value_line, plain_table, err)) {
                if (color) cout << CYAN;
                cout << "Results saved to " << outname << "\n" << (color?RESET:"");
                break;
            } else {
                if (color) cout << RED;
                cout << "Error writing file: " << err << "\n" << (color?RESET:"");
                cout << "Please enter a valid filename again.\n";
            }
        }
    }

    cout << "Done.\n";
    return 0;
}