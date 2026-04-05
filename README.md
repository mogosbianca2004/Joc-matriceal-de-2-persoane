# Joc-matriceal-de-2-persoane
#Teoria Jocurilor
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

using namespace std;

const double EPSILON = 1e-8;
const double BIG_M = 1000000.0;

// ================= SIMPLEX =================
class UniversalSimplex {
private:
    vector<vector<double>> tableau;
    vector<int> basic_variables;
    int num_rows, num_cols;
    int num_orig_vars;
    int first_artificial_idx;
    int num_artificial;
    bool is_max;

public:
    UniversalSimplex(vector<vector<double>> A, vector<double> b, vector<double> c, vector<int> signs, bool is_maximize) {
        is_max = is_maximize;
        int m = A.size();
        int n = c.size();
        num_orig_vars = n;

        for (int i = 0; i < m; ++i)
            if (b[i] < 0) {
                b[i] = -b[i];
                for (int j = 0; j < n; ++j) A[i][j] = -A[i][j];
                signs[i] = -signs[i];
            }

        int num_slack_surplus = 0;
        num_artificial = 0;

        for (int s : signs) {
            if (s != 0) num_slack_surplus++;
            if (s >= 0) num_artificial++;
        }

        first_artificial_idx = n + num_slack_surplus;
        num_rows = m + 1;
        num_cols = n + num_slack_surplus + num_artificial + 1;

        tableau.assign(num_rows, vector<double>(num_cols, 0.0));
        basic_variables.resize(m);

        int current_slack = n;
        int current_artificial = first_artificial_idx;

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j)
                tableau[i][j] = A[i][j];

            tableau[i][num_cols - 1] = b[i];

            if (signs[i] == -1) {
                tableau[i][current_slack] = 1.0;
                basic_variables[i] = current_slack++;
            }
            else if (signs[i] == 1) {
                tableau[i][current_slack] = -1.0;
                current_slack++;
                tableau[i][current_artificial] = 1.0;
                basic_variables[i] = current_artificial++;
            }
            else {
                tableau[i][current_artificial] = 1.0;
                basic_variables[i] = current_artificial++;
            }
        }

        for (int j = 0; j < n; ++j)
            tableau[m][j] = is_max ? -c[j] : c[j];

        int temp_artif = first_artificial_idx;
        for (int i = 0; i < m; ++i) {
            if (signs[i] >= 0) {
                tableau[m][temp_artif] = BIG_M;
                for (int j = 0; j < num_cols; ++j)
                    tableau[m][j] -= BIG_M * tableau[i][j];
                temp_artif++;
            }
        }
    }

    void solve(vector<double>& solution, double& optimal_value) {
        while (true) {
            int pivot_col = -1;
            double min_val = -EPSILON;

            for (int j = 0; j < num_cols - 1; ++j)
                if (tableau[num_rows - 1][j] < min_val) {
                    min_val = tableau[num_rows - 1][j];
                    pivot_col = j;
                }

            if (pivot_col == -1) break;

            int pivot_row = -1;
            double min_ratio = 1e30;

            for (int i = 0; i < num_rows - 1; ++i)
                if (tableau[i][pivot_col] > EPSILON) {
                    double ratio = tableau[i][num_cols - 1] / tableau[i][pivot_col];
                    if (ratio < min_ratio) {
                        min_ratio = ratio;
                        pivot_row = i;
                    }
                }

            if (pivot_row == -1)
                throw runtime_error("Solutie nemarginita!");

            double pivot_val = tableau[pivot_row][pivot_col];
            for (int j = 0; j < num_cols; ++j)
                tableau[pivot_row][j] /= pivot_val;

            for (int i = 0; i < num_rows; ++i)
                if (i != pivot_row) {
                    double factor = tableau[i][pivot_col];
                    for (int j = 0; j < num_cols; ++j)
                        tableau[i][j] -= factor * tableau[pivot_row][j];
                }
            basic_variables[pivot_row] = pivot_col;
        }

        solution.assign(num_orig_vars, 0.0);
        for (int i = 0; i < num_rows - 1; ++i)
            if (basic_variables[i] < num_orig_vars)
                solution[basic_variables[i]] = tableau[i][num_cols - 1];

        optimal_value = is_max ? tableau[num_rows - 1][num_cols - 1]
                               : -tableau[num_rows - 1][num_cols - 1];
    }
};

// ================= GAME THEORY =================

int main() {
    int m, n;
    cout << "Strategii A (linii): "; cin >> m;
    cout << "Strategii B (coloane): "; cin >> n;

    vector<vector<double>> Q(m, vector<double>(n));
    cout << "\nMatricea jocului:\n";
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            cout << "Q[" << i+1 << "][" << j+1 << "] = ";
            cin >> Q[i][j];
        }

    vector<double> alpha(m), beta(n);
    for (int i = 0; i < m; i++)
        alpha[i] = *min_element(Q[i].begin(), Q[i].end());

    for (int j = 0; j < n; j++) {
        double mx = Q[0][j];
        for (int i = 1; i < m; i++)
            mx = max(mx, Q[i][j]);
        beta[j] = mx;
    }

    double v1 = *max_element(alpha.begin(), alpha.end());
    double v2 = *min_element(beta.begin(), beta.end());

    // Setez afisarea la 2 zecimale
    cout << fixed << setprecision(2);
    
    cout << "\nAlpha: "; for (double x : alpha) cout << x << " ";
    cout << "\nBeta: "; for (double x : beta) cout << x << " ";
    cout << "\nv1 = " << v1 << " | v2 = " << v2 << endl;

    if (abs(v1 - v2) < EPSILON) {
        cout << "\n>>> PUNCT SA!\n";
        int linie = max_element(alpha.begin(), alpha.end()) - alpha.begin();
        int coloana = min_element(beta.begin(), beta.end()) - beta.begin();
        cout << "A joaca: a" << linie+1 << " | B joaca: b" << coloana+1 << endl;
        cout << "Valoare joc v = " << v1 << endl;
        return 0;
    }

    cout << "\n>>> STRATEGII MIXTE (Simplex)\n";

    // --- CALCUL STRATEGIA A ---
    vector<vector<double>> A_lp;
    vector<double> b_lp;
    vector<int> signs;

    for (int j = 0; j < n; j++) {
        vector<double> row(m + 1);
        for (int i = 0; i < m; i++) row[i] = Q[i][j];
        row[m] = -1; 
        A_lp.push_back(row);
        b_lp.push_back(0);
        signs.push_back(1); 
    }
    vector<double> sum_row_A(m + 1, 0);
    for (int i = 0; i < m; i++) sum_row_A[i] = 1;
    A_lp.push_back(sum_row_A); b_lp.push_back(1); signs.push_back(0); 

    vector<double> c_A(m + 1, 0); c_A[m] = 1; 

    // --- CALCUL STRATEGIA B ---
    vector<vector<double>> B_lp;
    vector<double> b_lp_B;
    vector<int> signs_B;

    for (int i = 0; i < m; i++) {
        vector<double> row(n + 1);
        for (int j = 0; j < n; j++) row[j] = Q[i][j];
        row[n] = -1; 
        B_lp.push_back(row);
        b_lp_B.push_back(0);
        signs_B.push_back(-1); 
    }
    vector<double> sum_row_B(n + 1, 0);
    for (int j = 0; j < n; j++) sum_row_B[j] = 1;
    B_lp.push_back(sum_row_B); b_lp_B.push_back(1); signs_B.push_back(0); 

    vector<double> c_B(n + 1, 0); c_B[n] = 1; 

    try {
        vector<double> solA, solB;
        double valA, valB;

        UniversalSimplex solverA(A_lp, b_lp, c_A, signs, true);
        solverA.solve(solA, valA);

        UniversalSimplex solverB(B_lp, b_lp_B, c_B, signs_B, false);
        solverB.solve(solB, valB);

        cout << "\nStrategia A (P): ";
        double sumA = 0;
        for (int i = 0; i < m; i++) {
            cout << "x" << i+1 << "=" << solA[i] << " ";
            sumA += solA[i];
        }
        cout << "| Suma A = " << sumA << endl;

        cout << "Strategia B (q): ";
        double sumB = 0;
        for (int j = 0; j < n; j++) {
            cout << "y" << j+1 << "=" << solB[j] << " ";
            sumB += solB[j];
        }
        cout << "| Suma B = " << sumB << endl;

        double v_final = solA[m];
        cout << "Valoare joc v = " << v_final << endl;

        // ===== VERIFICARE: v = P^T * Q * q =====
        double v_verif = 0;
        for (int i = 0; i < m; i++) {
            double temp = 0;
            for (int j = 0; j < n; j++) {
                temp += Q[i][j] * solB[j];
            }
            v_verif += solA[i] * temp;
        }

        cout << "\n>>> VERIFICARE FORMULA v = xo * Q * yo\n";
        cout << "v calculat prin formula: " << v_verif << endl;
 

    } catch (exception& e) {
        cout << "Eroare: " << e.what() << endl;
    }

    return 0;
}
