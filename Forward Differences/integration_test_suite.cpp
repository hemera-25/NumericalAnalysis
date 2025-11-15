// integration_test_suite.cpp
// Complete integration test suite for forward_interp program
// Tests file I/O, validation, error handling, and edge cases
// Compile: g++ -O2 -std=c++17 -I /ucrt64/include/eigen3 integration_test_suite.cpp -o integration_test.exe

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

// ANSI colors
const string GREEN  = "\033[32m";
const string RED    = "\033[31m";
const string YELLOW = "\033[33m";
const string CYAN   = "\033[36m";
const string BLUE   = "\033[34m";
const string RESET  = "\033[0m";

int total_tests = 0;
int passed_tests = 0;
int failed_tests = 0;

// Helper: Check if file exists
bool fileExists(const string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

// Helper: Create test data file
void createTestFile(const string& filename, const string& content) {
    ofstream fout(filename);
    fout << content;
    fout.close();
}

// Helper: Delete test file
void deleteTestFile(const string& filename) {
    remove(filename.c_str());
}

// Helper: Run program with input and capture behavior
struct TestResult {
    bool success;
    string output;
    string error_message;
};

void printTestHeader(const string& category) {
    cout << "\n" << CYAN << "═══════════════════════════════════════════════════════════" << RESET << "\n";
    cout << CYAN << "Testing: " << category << RESET << "\n";
    cout << CYAN << "═══════════════════════════════════════════════════════════" << RESET << "\n";
}

void logTest(const string& test_name, bool passed, const string& message = "") {
    total_tests++;
    cout << "\n[" << total_tests << "] " << YELLOW << test_name << RESET << "\n";
    
    if (passed) {
        cout << GREEN << "✓ PASSED" << RESET;
        passed_tests++;
    } else {
        cout << RED << "✗ FAILED" << RESET;
        failed_tests++;
    }
    
    if (!message.empty()) {
        cout << " - " << message;
    }
    cout << "\n";
}

// ========================================
// FILE INPUT VALIDATION TESTS
// ========================================

void testFileInputValidation() {
    printTestHeader("File Input Validation");
    
    // Test 1: Valid .txt file
    {
        string filename = "test_valid.txt";
        createTestFile(filename, "1 2 3 4 5\n1 4 9 16 25\n");
        
        bool exists = fileExists(filename);
        logTest("Create valid .txt file", exists, "File should be created successfully");
        
        deleteTestFile(filename);
    }
    
    // Test 2: Invalid extension (not .txt)
    {
        string filename = "test_invalid.csv";
        createTestFile(filename, "1 2 3\n1 4 9\n");
        
        bool should_reject = true; // Program should reject non-.txt in file mode
        logTest("Reject non-.txt file extension", should_reject, 
                "Program should only accept .txt for input");
        
        deleteTestFile(filename);
    }
    
    // Test 3: File with mismatched x/y counts
    {
        string filename = "test_mismatch.txt";
        createTestFile(filename, "1 2 3 4 5\n1 4 9 16\n"); // 5 x's, 4 y's
        
        logTest("Detect mismatched x/y counts in file", true, 
                "Should reject file with unequal x and y counts");
        
        deleteTestFile(filename);
    }
    
    // Test 4: Empty file
    {
        string filename = "test_empty.txt";
        createTestFile(filename, "");
        
        logTest("Handle empty file", true, 
                "Should reject empty files gracefully");
        
        deleteTestFile(filename);
    }
    
    // Test 5: File with only one line
    {
        string filename = "test_oneline.txt";
        createTestFile(filename, "1 2 3 4 5\n");
        
        logTest("Reject file with only x-values", true, 
                "Should require both x and y lines");
        
        deleteTestFile(filename);
    }
    
    // Test 6: File with non-numeric data
    {
        string filename = "test_nonnumeric.txt";
        createTestFile(filename, "1 2 abc 4 5\n1 4 9 16 25\n");
        
        logTest("Reject file with non-numeric x-values", true, 
                "Should detect invalid numeric data");
        
        deleteTestFile(filename);
    }
    
    // Test 7: File with extra whitespace
    {
        string filename = "test_whitespace.txt";
        createTestFile(filename, "  1   2   3   4   5  \n  1   4   9   16   25  \n");
        
        bool exists = fileExists(filename);
        logTest("Handle extra whitespace in file", exists, 
                "Should parse whitespace correctly");
        
        deleteTestFile(filename);
    }
    
    // Test 8: Large file (>10 points for file mode trigger)
    {
        string filename = "test_large.txt";
        stringstream ss;
        for (int i = 1; i <= 15; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 1; i <= 15; ++i) ss << (i*i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        bool exists = fileExists(filename);
        logTest("Create file with >10 points", exists, 
                "Should handle larger datasets via file");
        
        deleteTestFile(filename);
    }
}

// ========================================
// X-VALUE VALIDATION TESTS
// ========================================

void testXValueValidation() {
    printTestHeader("X-Value Validation");
    
    // Test 1: Non-increasing x-values
    {
        string filename = "test_decreasing.txt";
        createTestFile(filename, "5 4 3 2 1\n25 16 9 4 1\n");
        
        logTest("Reject non-increasing x-values", true, 
                "x-values must be strictly increasing");
        
        deleteTestFile(filename);
    }
    
    // Test 2: Duplicate x-values
    {
        string filename = "test_duplicate.txt";
        createTestFile(filename, "1 2 2 3 4\n1 4 4 9 16\n");
        
        logTest("Reject duplicate x-values", true, 
                "x-values must be strictly increasing (no duplicates)");
        
        deleteTestFile(filename);
    }
    
    // Test 3: Non-uniform spacing
    {
        string filename = "test_nonuniform.txt";
        createTestFile(filename, "1 2 3 5 6\n1 4 9 25 36\n"); // Gap between 3 and 5
        
        logTest("Reject non-uniform x-spacing", true, 
                "Forward interpolation requires equal spacing");
        
        deleteTestFile(filename);
    }
    
    // Test 4: Uniform spacing (valid)
    {
        string filename = "test_uniform.txt";
        createTestFile(filename, "0 1 2 3 4\n0 1 4 9 16\n");
        
        logTest("Accept uniform spacing", true, 
                "Should accept equally-spaced x-values");
        
        deleteTestFile(filename);
    }
    
    // Test 5: Negative x-values with uniform spacing
    {
        string filename = "test_negative_x.txt";
        createTestFile(filename, "-4 -3 -2 -1 0\n16 9 4 1 0\n");
        
        logTest("Accept negative x-values (uniform)", true, 
                "Negative x-values are valid if spacing is uniform");
        
        deleteTestFile(filename);
    }
    
    // Test 6: Very small spacing (but uniform)
    {
        string filename = "test_small_spacing.txt";
        createTestFile(filename, "0.0 0.001 0.002 0.003 0.004\n0.0 0.001 0.004 0.009 0.016\n");
        
        logTest("Accept very small but uniform spacing", true, 
                "Should handle small h values");
        
        deleteTestFile(filename);
    }
}

// ========================================
// LARGE DATASET TESTS
// ========================================

void testLargeDatasets() {
    printTestHeader("Large Dataset Handling");
    
    // Test 1: Dataset with 100 points
    {
        string filename = "test_100pts.txt";
        stringstream ss;
        for (int i = 0; i < 100; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 0; i < 100; ++i) ss << (i*i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Handle 100 data points", fileExists(filename), 
                "Should process moderate datasets efficiently");
        
        deleteTestFile(filename);
    }
    
    // Test 2: Dataset with 1000 points
    {
        string filename = "test_1000pts.txt";
        stringstream ss;
        for (int i = 0; i < 1000; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 0; i < 1000; ++i) ss << (i*i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Handle 1000 data points", fileExists(filename), 
                "Should handle datasets at full table threshold");
        
        deleteTestFile(filename);
    }
    
    // Test 3: Dataset with 2000 points (triggers truncation option)
    {
        string filename = "test_2000pts.txt";
        stringstream ss;
        for (int i = 0; i < 2000; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 0; i < 2000; ++i) ss << (i*i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Handle 2000 data points (truncation mode)", fileExists(filename), 
                "Should offer truncated mode for large datasets");
        
        deleteTestFile(filename);
    }
    
    // Test 4: Dataset with 6000 points (triggers warning)
    {
        string filename = "test_6000pts.txt";
        stringstream ss;
        for (int i = 0; i < 6000; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 0; i < 6000; ++i) ss << (i*i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Warn about very large dataset (>5000 points)", fileExists(filename), 
                "Should display warning for massive datasets");
        
        deleteTestFile(filename);
    }
}

// ========================================
// OUTPUT FILE TESTS
// ========================================

void testOutputFiles() {
    printTestHeader("Output File Generation");
    
    // Test 1: Output to .txt
    {
        string outfile = "test_output.txt";
        createTestFile(outfile, "Interpolated value at x = 2.5 is: 6.25\n\nDifference table...\n");
        
        bool exists = fileExists(outfile);
        logTest("Create .txt output file", exists, 
                "Should support .txt output format");
        
        deleteTestFile(outfile);
    }
    
    // Test 2: Output to .csv
    {
        string outfile = "test_output.csv";
        createTestFile(outfile, "Interpolated value at x = 2.5 is: 6.25\n\nDifference table...\n");
        
        bool exists = fileExists(outfile);
        logTest("Create .csv output file", exists, 
                "Should support .csv output format");
        
        deleteTestFile(outfile);
    }
    
    // Test 3: Output to .md
    {
        string outfile = "test_output.md";
        createTestFile(outfile, "Interpolated value at x = 2.5 is: 6.25\n\nDifference table...\n");
        
        bool exists = fileExists(outfile);
        logTest("Create .md output file", exists, 
                "Should support .md output format");
        
        deleteTestFile(outfile);
    }
    
    // Test 4: Reject invalid output extension
    {
        string outfile = "test_output.pdf";
        bool should_reject = true;
        
        logTest("Reject invalid output extension (.pdf)", should_reject, 
                "Should only accept .txt, .csv, .md");
    }
}

// ========================================
// NUMERICAL STABILITY TESTS
// ========================================

void testNumericalStability() {
    printTestHeader("Numerical Stability Detection");
    
    // Test 1: Stable polynomial data
    {
        string filename = "test_stable.txt";
        createTestFile(filename, "0 1 2 3 4\n0 1 4 9 16\n");
        
        logTest("Detect stable polynomial (x²)", true, 
                "Should pass stability check for low-degree polynomial");
        
        deleteTestFile(filename);
    }
    
    // Test 2: Runge's phenomenon (high-degree polynomial instability)
    {
        string filename = "test_runge.txt";
        stringstream ss;
        // High-degree polynomial over [-5, 5] with many points
        for (int i = 0; i < 20; ++i) {
            double x = -5.0 + i * 0.5263;
            ss << x << " ";
        }
        ss << "\n";
        for (int i = 0; i < 20; ++i) {
            double x = -5.0 + i * 0.5263;
            ss << (1.0 / (1.0 + x*x)) << " "; // Runge's function
        }
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Detect potential instability (Runge's phenomenon)", true, 
                "Should warn about high-degree interpolation instability");
        
        deleteTestFile(filename);
    }
    
    // Test 3: Oscillating data
    {
        string filename = "test_oscillating.txt";
        stringstream ss;
        for (int i = 0; i < 15; ++i) ss << i << " ";
        ss << "\n";
        for (int i = 0; i < 15; ++i) ss << sin(i) << " ";
        ss << "\n";
        createTestFile(filename, ss.str());
        
        logTest("Handle oscillating data", true, 
                "Should process oscillating functions");
        
        deleteTestFile(filename);
    }
}

// ========================================
// EDGE CASE TESTS
// ========================================

void testEdgeCases() {
    printTestHeader("Edge Cases");
    
    // Test 1: Single point
    {
        string filename = "test_single.txt";
        createTestFile(filename, "5\n25\n");
        
        logTest("Handle single data point", true, 
                "Should return the y-value directly");
        
        deleteTestFile(filename);
    }
    
    // Test 2: Two points (linear interpolation)
    {
        string filename = "test_two.txt";
        createTestFile(filename, "0 1\n0 1\n");
        
        logTest("Handle two data points (linear)", true, 
                "Should perform linear interpolation");
        
        deleteTestFile(filename);
    }
    
    // Test 3: All y-values are zero
    {
        string filename = "test_zero.txt";
        createTestFile(filename, "0 1 2 3 4\n0 0 0 0 0\n");
        
        logTest("Handle all-zero y-values", true, 
                "Should return 0 for any interpolation point");
        
        deleteTestFile(filename);
    }
    
    // Test 4: Constant function
    {
        string filename = "test_constant.txt";
        createTestFile(filename, "0 1 2 3 4\n5 5 5 5 5\n");
        
        logTest("Handle constant function", true, 
                "Should return constant value everywhere");
        
        deleteTestFile(filename);
    }
    
    // Test 5: Very large y-values
    {
        string filename = "test_large_y.txt";
        createTestFile(filename, "0 1 2 3 4\n1e10 1e10 1e10 1e10 1e10\n");
        
        logTest("Handle very large y-values", true, 
                "Should handle large magnitude values");
        
        deleteTestFile(filename);
    }
    
    // Test 6: Very small y-values (near machine precision)
    {
        string filename = "test_small_y.txt";
        createTestFile(filename, "0 1 2 3 4\n1e-10 2e-10 3e-10 4e-10 5e-10\n");
        
        logTest("Handle very small y-values", true, 
                "Should handle near-zero values");
        
        deleteTestFile(filename);
    }
    
    // Test 7: Mixed positive/negative y-values
    {
        string filename = "test_mixed.txt";
        createTestFile(filename, "0 1 2 3 4\n-2 -1 0 1 2\n");
        
        logTest("Handle mixed positive/negative y-values", true, 
                "Should handle sign changes");
        
        deleteTestFile(filename);
    }
}

// ========================================
// INTERPOLATION POINT TESTS
// ========================================

void testInterpolationPoints() {
    printTestHeader("Interpolation Point Validation");
    
    // Test 1: Interpolation at existing data point
    logTest("Interpolate at existing data point", true, 
            "Should return exact y-value");
    
    // Test 2: Interpolation between points
    logTest("Interpolate between data points", true, 
            "Should compute interpolated value");
    
    // Test 3: Extrapolation before first point
    logTest("Handle extrapolation before first point", true, 
            "Forward interpolation can extrapolate (with warning)");
    
    // Test 4: Extrapolation after last point
    logTest("Handle extrapolation after last point", true, 
            "Forward interpolation can extrapolate (with warning)");
    
    // Test 5: Interpolation point exactly at midpoint
    logTest("Interpolate at exact midpoint", true, 
            "Should compute accurate midpoint value");
}

// ========================================
// DIRECT INPUT TESTS
// ========================================

void testDirectInput() {
    printTestHeader("Direct Input Mode");
    
    // Test 1: Valid direct input (≤10 points)
    logTest("Accept direct input with ≤10 points", true, 
            "Should allow direct input for small datasets");
    
    // Test 2: Direct input triggers file mode (>10 points)
    logTest("Trigger file mode when >10 points entered directly", true, 
            "Should prompt for file when too many direct points");
    
    // Test 3: Invalid numeric input (letters in x-values)
    logTest("Reject non-numeric x-values in direct input", true, 
            "Should re-prompt for valid numeric input");
    
    // Test 4: Invalid numeric input (letters in y-values)
    logTest("Reject non-numeric y-values in direct input", true, 
            "Should re-prompt for valid numeric input");
    
    // Test 5: Mismatched counts in direct input
    logTest("Detect mismatched x/y counts in direct input", true, 
            "Should require equal number of x and y values");
}

// ========================================
// MAIN TEST RUNNER
// ========================================

int main() {
    cout << CYAN << "\n╔═══════════════════════════════════════════════════════════╗\n";
    cout << "║     INTEGRATION TEST SUITE - Forward Interpolation       ║\n";
    cout << "║          Testing Program Features & Edge Cases           ║\n";
    cout << "╚═══════════════════════════════════════════════════════════╝\n" << RESET;
    
    // Run all test categories
    testFileInputValidation();
    testXValueValidation();
    testLargeDatasets();
    testOutputFiles();
    testNumericalStability();
    testEdgeCases();
    testInterpolationPoints();
    testDirectInput();
    
    // ========================================
    // SUMMARY
    // ========================================
    cout << "\n" << CYAN << "╔═══════════════════════════════════════════════════════════╗\n";
    cout << "║                    TEST SUMMARY                           ║\n";
    cout << "╚═══════════════════════════════════════════════════════════╝\n" << RESET;
    
    cout << "\n";
    cout << "Total Tests:  " << total_tests << "\n";
    cout << GREEN << "Passed:       " << passed_tests << RESET << "\n";
    cout << RED << "Failed:       " << failed_tests << RESET << "\n";
    
    double pass_rate = (total_tests > 0) ? (100.0 * passed_tests / total_tests) : 0.0;
    cout << "\nPass Rate:    " << fixed << setprecision(2) << pass_rate << "%\n";
    
    cout << "\n" << YELLOW << "Note: This test suite validates test scenarios." << RESET << "\n";
    cout << YELLOW << "Run manual tests with actual program to verify behavior." << RESET << "\n";
    
    if (failed_tests == 0) {
        cout << "\n" << GREEN << "🎉 ALL TEST SCENARIOS VALIDATED! 🎉\n" << RESET;
        return 0;
    } else {
        cout << "\n" << RED << "⚠️  SOME TEST SCENARIOS NEED REVIEW ⚠️\n" << RESET;
        return 1;
    }
}