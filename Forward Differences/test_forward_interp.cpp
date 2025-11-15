// test_forward_interp.cpp
// Comprehensive automated test suite for Newton-Gregory Forward Interpolation
// Compile: g++ -O2 -std=c++17 -I /ucrt64/include/eigen3 test_forward_interp.cpp -o test_suite.exe

#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <functional>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ANSI color codes
const string GREEN  = "\033[32m";
const string RED    = "\033[31m";
const string YELLOW = "\033[33m";
const string CYAN   = "\033[36m";
const string RESET  = "\033[0m";

// Test statistics
int total_tests = 0;
int passed_tests = 0;
int failed_tests = 0;

// Forward interpolation implementation (copied from main program)
double factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) result *= i;
    return result;
}

MatrixXd buildForwardTable(const vector<double> &y) {
    int n = (int)y.size();
    MatrixXd D = MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) D(i,0) = y[i];
    for (int j = 1; j < n; ++j)
        for (int i = 0; i + j < n; ++i)
            D(i,j) = D(i+1,j-1) - D(i,j-1);
    return D;
}

double interpolateNewtonForward(const MatrixXd &D, const vector<double> &x, double xp) {
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

// Test framework
struct TestCase {
    string name;
    vector<double> x;
    vector<double> y;
    double xp;
    double expected;
    double tolerance;
    bool should_pass;
    string category;
};

void printTestHeader(const string& category) {
    cout << "\n" << CYAN << "========================================" << RESET << "\n";
    cout << CYAN << "Testing: " << category << RESET << "\n";
    cout << CYAN << "========================================" << RESET << "\n";
}

void runTest(const TestCase& test) {
    total_tests++;
    cout << "\n[" << total_tests << "] " << YELLOW << test.name << RESET << "\n";
    
    try {
        // Build difference table
        MatrixXd D = buildForwardTable(test.y);
        
        // Perform interpolation
        double result = interpolateNewtonForward(D, test.x, test.xp);
        
        // Check result
        double error = fabs(result - test.expected);
        bool passed = error <= test.tolerance;
        
        if (passed && test.should_pass) {
            cout << GREEN << "✓ PASSED" << RESET << "\n";
            cout << "  Expected: " << fixed << setprecision(6) << test.expected << "\n";
            cout << "  Got:      " << result << "\n";
            cout << "  Error:    " << scientific << error << "\n";
            passed_tests++;
        } else if (!passed && !test.should_pass) {
            cout << GREEN << "✓ PASSED (expected failure)" << RESET << "\n";
            passed_tests++;
        } else {
            cout << RED << "✗ FAILED" << RESET << "\n";
            cout << "  Expected: " << fixed << setprecision(6) << test.expected << "\n";
            cout << "  Got:      " << result << "\n";
            cout << "  Error:    " << scientific << error << " (tolerance: " << test.tolerance << ")\n";
            failed_tests++;
        }
    } catch (const exception& e) {
        if (!test.should_pass) {
            cout << GREEN << "✓ PASSED (caught expected exception)" << RESET << "\n";
            passed_tests++;
        } else {
            cout << RED << "✗ FAILED (exception)" << RESET << "\n";
            cout << "  Exception: " << e.what() << "\n";
            failed_tests++;
        }
    }
}

// Generate test data from function
vector<TestCase> generateFunctionTests(const string& func_name, 
                                       function<double(double)> f,
                                       double x_start, double x_end, int n_points,
                                       double tolerance) {
    vector<TestCase> tests;
    
    // Generate equally-spaced points
    double h = (x_end - x_start) / (n_points - 1);
    vector<double> x(n_points), y(n_points);
    
    for (int i = 0; i < n_points; ++i) {
        x[i] = x_start + i * h;
        y[i] = f(x[i]);
    }
    
    // Test interpolation at various points
    vector<pair<string, double>> test_points = {
        {"at data point", x[n_points/2]},
        {"midpoint between nodes", x[0] + h/2},
        {"quarter point", x[0] + h/4},
        {"near start", x[0] + h/10},
        {"near end", x[n_points-1] - h/10}
    };
    
    for (const auto& [desc, xp] : test_points) {
        TestCase test;
        test.name = func_name + " - " + desc;
        test.x = x;
        test.y = y;
        test.xp = xp;
        test.expected = f(xp);
        test.tolerance = tolerance;
        test.should_pass = true;
        test.category = func_name;
        tests.push_back(test);
    }
    
    return tests;
}

int main() {
    cout << CYAN << "\n╔═══════════════════════════════════════════════════════════╗\n";
    cout << "║  Newton-Gregory Forward Interpolation Test Suite         ║\n";
    cout << "╚═══════════════════════════════════════════════════════════╝\n" << RESET;
    
    vector<TestCase> all_tests;
    
    // ========================================
    // 1. EDGE CASES
    // ========================================
    printTestHeader("Edge Cases");
    
    // Single point
    all_tests.push_back({
        "Single data point",
        {5.0}, {25.0}, 5.0, 25.0, 1e-10, true, "Edge Cases"
    });
    
    // Two points (linear)
    all_tests.push_back({
        "Two points - exact interpolation",
        {0.0, 1.0}, {0.0, 1.0}, 0.5, 0.5, 1e-10, true, "Edge Cases"
    });
    
    // Three points (parabola exact fit)
    all_tests.push_back({
        "Three points - parabola x²",
        {0.0, 1.0, 2.0}, {0.0, 1.0, 4.0}, 1.5, 2.25, 1e-10, true, "Edge Cases"
    });
    
    // ========================================
    // 2. POLYNOMIAL FUNCTIONS (Exact Fit)
    // ========================================
    printTestHeader("Polynomial Functions - Exact Fit");
    
    // Linear: y = 2x + 3
    auto linear_tests = generateFunctionTests(
        "Linear (y=2x+3)", 
        [](double x) { return 2*x + 3; },
        0.0, 10.0, 5, 1e-10
    );
    all_tests.insert(all_tests.end(), linear_tests.begin(), linear_tests.end());
    
    // Quadratic: y = x²
    auto quad_tests = generateFunctionTests(
        "Quadratic (y=x²)",
        [](double x) { return x*x; },
        0.0, 5.0, 6, 1e-8
    );
    all_tests.insert(all_tests.end(), quad_tests.begin(), quad_tests.end());
    
    // Cubic: y = x³ - 2x² + 3x - 1
    auto cubic_tests = generateFunctionTests(
        "Cubic (y=x³-2x²+3x-1)",
        [](double x) { return x*x*x - 2*x*x + 3*x - 1; },
        0.0, 4.0, 8, 1e-6
    );
    all_tests.insert(all_tests.end(), cubic_tests.begin(), cubic_tests.end());
    
    // Quartic: y = x⁴
    auto quartic_tests = generateFunctionTests(
        "Quartic (y=x⁴)",
        [](double x) { return x*x*x*x; },
        0.0, 3.0, 10, 1e-4
    );
    all_tests.insert(all_tests.end(), quartic_tests.begin(), quartic_tests.end());
    
    // ========================================
    // 3. TRANSCENDENTAL FUNCTIONS (Approximation)
    // ========================================
    printTestHeader("Transcendental Functions - Approximation");
    
    // Exponential: y = e^x
    auto exp_tests = generateFunctionTests(
        "Exponential (y=e^x)",
        [](double x) { return exp(x); },
        0.0, 2.0, 10, 0.01
    );
    all_tests.insert(all_tests.end(), exp_tests.begin(), exp_tests.end());
    
    // Sine: y = sin(x)
    auto sin_tests = generateFunctionTests(
        "Sine (y=sin(x))",
        [](double x) { return sin(x); },
        0.0, M_PI, 12, 0.001
    );
    all_tests.insert(all_tests.end(), sin_tests.begin(), sin_tests.end());
    
    // Cosine: y = cos(x)
    auto cos_tests = generateFunctionTests(
        "Cosine (y=cos(x))",
        [](double x) { return cos(x); },
        0.0, M_PI, 12, 0.001
    );
    all_tests.insert(all_tests.end(), cos_tests.begin(), cos_tests.end());
    
    // Logarithm: y = ln(x+1)
    auto log_tests = generateFunctionTests(
        "Logarithm (y=ln(x+1))",
        [](double x) { return log(x + 1); },
        0.0, 5.0, 15, 0.01
    );
    all_tests.insert(all_tests.end(), log_tests.begin(), log_tests.end());
    
    // ========================================
    // 4. SPECIAL FUNCTIONS
    // ========================================
    printTestHeader("Special Functions");
    
    // Square root: y = √x
    auto sqrt_tests = generateFunctionTests(
        "Square Root (y=√x)",
        [](double x) { return sqrt(x); },
        1.0, 10.0, 15, 0.01
    );
    all_tests.insert(all_tests.end(), sqrt_tests.begin(), sqrt_tests.end());
    
    // Reciprocal: y = 1/x
    auto reciprocal_tests = generateFunctionTests(
        "Reciprocal (y=1/x)",
        [](double x) { return 1.0/x; },
        1.0, 5.0, 12, 0.01
    );
    all_tests.insert(all_tests.end(), reciprocal_tests.begin(), reciprocal_tests.end());
    
    // Gaussian: y = e^(-x²)
    auto gaussian_tests = generateFunctionTests(
        "Gaussian (y=e^(-x²))",
        [](double x) { return exp(-x*x); },
        -2.0, 2.0, 15, 0.01
    );
    all_tests.insert(all_tests.end(), gaussian_tests.begin(), gaussian_tests.end());
    
    // ========================================
    // 5. OSCILLATING FUNCTIONS
    // ========================================
    printTestHeader("Oscillating Functions");
    
    // Damped sine: y = e^(-x/2) * sin(2x)
    auto damped_sin_tests = generateFunctionTests(
        "Damped Sine (y=e^(-x/2)*sin(2x))",
        [](double x) { return exp(-x/2.0) * sin(2*x); },
        0.0, 2*M_PI, 20, 0.05
    );
    all_tests.insert(all_tests.end(), damped_sin_tests.begin(), damped_sin_tests.end());
    
    // ========================================
    // 6. CONSTANT FUNCTIONS
    // ========================================
    printTestHeader("Constant Function");
    
    all_tests.push_back({
        "Constant function (y=5)",
        {0.0, 1.0, 2.0, 3.0, 4.0},
        {5.0, 5.0, 5.0, 5.0, 5.0},
        2.5, 5.0, 1e-10, true, "Constant"
    });
    
    // ========================================
    // 7. NEGATIVE VALUES
    // ========================================
    printTestHeader("Negative Values");
    
    all_tests.push_back({
        "Negative x-values",
        {-4.0, -3.0, -2.0, -1.0, 0.0},
        {16.0, 9.0, 4.0, 1.0, 0.0},
        -2.5, 6.25, 1e-8, true, "Negative"
    });
    
    all_tests.push_back({
        "Negative y-values",
        {0.0, 1.0, 2.0, 3.0, 4.0},
        {0.0, -1.0, -4.0, -9.0, -16.0},
        2.5, -6.25, 1e-8, true, "Negative"
    });
    
    // ========================================
    // 8. LARGE SCALE VALUES
    // ========================================
    printTestHeader("Large Scale Values");
    
    all_tests.push_back({
        "Large x and y values",
        {1000.0, 1001.0, 1002.0, 1003.0, 1004.0},
        {1e6, 1.002001e6, 1.004004e6, 1.006009e6, 1.008016e6},
        1002.5, 1.00500625e6, 1.0, true, "Large Scale"
    });
    
    // ========================================
    // 9. SMALL SPACING
    // ========================================
    printTestHeader("Small Spacing");
    
    all_tests.push_back({
        "Very small spacing (h=0.001)",
        {0.0, 0.001, 0.002, 0.003, 0.004},
        {0.0, 0.001, 0.004, 0.009, 0.016},
        0.0025, 0.00625, 1e-8, true, "Small Spacing"
    });
    
    // ========================================
    // RUN ALL TESTS
    // ========================================
    
    string current_category = "";
    for (const auto& test : all_tests) {
        if (test.category != current_category) {
            current_category = test.category;
            if (current_category != "Edge Cases" && 
                current_category != "Polynomial Functions - Exact Fit" &&
                current_category != "Transcendental Functions - Approximation" &&
                current_category != "Special Functions" &&
                current_category != "Oscillating Functions" &&
                current_category != "Constant" &&
                current_category != "Negative" &&
                current_category != "Large Scale" &&
                current_category != "Small Spacing") {
                // Don't print header for function-generated tests
            }
        }
        runTest(test);
    }
    
    // ========================================
    // SUMMARY
    // ========================================
    cout << "\n" << CYAN << "╔═══════════════════════════════════════════════════════════╗\n";
    cout << "║                      TEST SUMMARY                         ║\n";
    cout << "╚═══════════════════════════════════════════════════════════╝\n" << RESET;
    
    cout << "\n";
    cout << "Total Tests:  " << total_tests << "\n";
    cout << GREEN << "Passed:       " << passed_tests << RESET << "\n";
    cout << RED << "Failed:       " << failed_tests << RESET << "\n";
    
    double pass_rate = (total_tests > 0) ? (100.0 * passed_tests / total_tests) : 0.0;
    cout << "\nPass Rate:    " << fixed << setprecision(2) << pass_rate << "%\n";
    
    if (failed_tests == 0) {
        cout << "\n" << GREEN << "🎉 ALL TESTS PASSED! 🎉\n" << RESET;
        return 0;
    } else {
        cout << "\n" << RED << "⚠️  SOME TESTS FAILED ⚠️\n" << RESET;
        return 1;
    }
}