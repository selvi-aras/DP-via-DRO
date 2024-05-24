// let me try my best
#include <iostream>        //input output
#include <vector>        //to add vector<type> variables
#include <chrono>        //to measure time ellapsed in a function
#include "gurobi_c++.h" //Gurobi's functions
#include <algorithm>
#include <unordered_map>
#include <math.h>       //some maths functions
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;
using namespace std::chrono;

GRBEnv* env; //GUROBI environment
GRBModel* model;//model
vector<vector<GRBVar>> p;
GRBConstr* c;

double epsilon;
double delta;

//new
std::vector<double> read_csv_to_vector(const std::string& filename) {
    std::vector<double> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    while (std::getline(file, line)) {
        std::stringstream linestream(line);
        double value;
        while (linestream >> value) {
            data.push_back(value);
            if (linestream.peek() == ',')
                linestream.ignore();
        }
    }

    file.close();
    return data;
}

/*
**************************************
PRELIMINARY Functions
**************************************
*/
void write_to_csv(const vector<double>& vect, ofstream& f) { // takes a double vector and makes it a CSV row of a given file
    f << vect[0];
    for (int i = 1; i < vect.size(); i++) {
        f << ',' << vect[i];
    }
    f << '\n';
}
void print_vector(const vector<double>& v) { //prints a <double> vector //cout << "Elements of the given vector are: " << endl;
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << '\n';
    }
}
void remove(vector<double>& v){//removes duplicates from a given vector -- O(nlog(n))
    sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
}

void sortVector(std::vector<double>& v) {
    std::sort(v.begin(), v.end());
}

void upper_lower(vector<double>& v, const double threshold) { //removes elements > upper_f || < - upper_f, and adds these two limits
    v.erase(std::remove_if(
        v.begin(), v.end(),
        [&](const double& x) {
            return x >= threshold || x <= -threshold; // condition s.t. we should remove from vector
        }), v.end());
    //now add end_points
    v.push_back(threshold);
    v.push_back(-threshold);
}
vector<double> differences(const vector<double>& v) {//{x_i - x_j \ : \ i,j = 1, ..., N+1}
    int size = (int)v.size();
    vector<double> diff; diff.reserve(size * size); //difference vector to return
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            diff.push_back( round((v[i] - v[j])*pow(10,5))/pow(10,5) );
        }
    }
    return diff;
}

void arange(vector<double>& result, double start, double end, double increase) { //Numpy's arange copied, but limits included, so small modification
    result.clear();
    //result.reserve((end - start) / increase);
    for (double i = start; i <= end + pow(10,-5); i = i + increase) {
        result.push_back(round(i * pow(10,7)) / pow(10,7));
    }
}

void add_dp_constr(GRBLinExpr& lhs, const vector<double>& X, int N, vector<vector<double>>& indexD, vector<vector<double>>& indexDD,vector<vector<double>>& indexD2, vector<vector<double>>& indexDD2, double& index, double& indexd, double& index2, double& indexd2) {
    lhs = 0.0; //start building the lhs of constraint
    for (int i = 0; i < 2 * N + 1; i++) {
        for (int j = 0; j < 2 * N + 1; j++){
            index = indexD[i][j];
            index2 = indexD2[i][j];
            
            indexd = indexDD[i][j];
            indexd2 = indexDD2[i][j];
            
//            cout << index << " ; " << index2 <<  " ; " << indexd <<  " ; " << indexd2 << endl;
            
            if (index >= 0 && index2 >= 0){
                lhs += p[index][index2];
            }
            if (indexd >= 0 && indexd2 >= 0){
                lhs -= exp(epsilon)*p[indexd][indexd2];
            }
            
        }
    }
    model->addConstr(lhs <= delta);
}
/*
**************************************
ESSENTIAL Functions Below (Algorithms)
**************************************
*/
void worst_event_improved(const double epsilon, const double varphi, const double varphi2, const vector<double>& X,  const vector<vector<double>>& p_decisions, double& violation){
    int N = (int)X.size() - 1; // size of p just to iterate over intervals
    
    double incr = round((X[1] - X[0])*pow(10,5))/pow(10,5); // take increments
    int shift_factor_1 = round(varphi/incr);
    int shift_factor_2 = round(varphi2/incr);
    int i_prime;
    int j_prime;
    violation = 0.0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            i_prime = i - shift_factor_1;
            j_prime = j - shift_factor_2;
            
            if (min(i_prime, j_prime) < 0 || max(i_prime, j_prime) >= N){
                violation += p_decisions[i][j];
            }
            else{
            
                violation += max(0.0, p_decisions[i][j] - exp(epsilon)*p_decisions[i_prime][j_prime]);
            }
        }
    }
}

void worst_event_to_add(const double epsilon, const double varphi, const double varphi2, const vector<double>& X,  const vector<vector<GRBVar>>& p, double& violation, GRBLinExpr& possible_worst_case){
    int N = (int)X.size() - 1; // size of p just to iterate over intervals
    
    violation = 0.0;
    possible_worst_case = 0;
    
    double incr = round((X[1] - X[0])*pow(10,5))/pow(10,5); // take increments
    int shift_factor_1 = round(varphi/incr);
    int shift_factor_2 = round(varphi2/incr);
    int i_prime;
    int j_prime;
    
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            i_prime = i - shift_factor_1;
            j_prime = j - shift_factor_2;
            
            if (min(i_prime, j_prime) < 0 || max(i_prime, j_prime) >= N){
                violation += p[i][j].get(GRB_DoubleAttr_X);
                possible_worst_case += p[i][j];
            }
            else{
                if (p[i][j].get(GRB_DoubleAttr_X) - exp(epsilon)*p[i_prime][j_prime].get(GRB_DoubleAttr_X) >= pow(10,-8)){
                    violation += p[i][j].get(GRB_DoubleAttr_X) - exp(epsilon)*p[i_prime][j_prime].get(GRB_DoubleAttr_X);
                    possible_worst_case += p[i][j] - exp(epsilon)*p[i_prime][j_prime];
                }
                
            }
        }
    }
}

//ALG3 below
void worst_event(const double epsilon, const double varphi,const double varphi2, const vector<double>& X, const vector<vector<GRBVar>>& p, vector<vector<double>>& indexD, vector<vector<double>>& indexD2, vector<vector<double>>& indexDD, vector<vector<double>>& indexDD2, double& violation) {
    //returns the worst-case violation for a fixed varphi (ALG3 of the paper!)
    //assume wlog that f(D) = 0 and f(D') = varphi
    
    //start with the first dimension now
    double condition;//condition to check
    double length_ell; double length_ell2; //length of each interval
    bool check_if_i; bool check_if_i2; //checks if the interval starts with an f(D)+i point or not
    violation = 0; //new addition (normally 0) //initially no violation
    int N = (int)X.size() - 1; // -1 because X has size N+1 and probability vector has size N (N intervals)
    int i = -1, j = -1; //position of the last noise index f(D) + x_i or f(D') + x_j. =-1 means no access or don't include this interval
    int i2, j2;
    double f_under_ell; //this will keep the beginning of A_ell = [f_under_ell, f_over_ell], one interval from the common refinement
    double f_under_ell2; //this will keep the beginning of A_ell = [f_under_ell, f_over_ell], one interval from the common refinement
    //****** Initialize the first f_underline_ell
    if (varphi > 0) {//means f(D) < f(D') so f(D) + x_i comes first in the common refinement
        i++;
        f_under_ell = X[i];
    }
    else {
        j++;
        f_under_ell = X[j] + varphi;
    }
    
    
    //****** main for loop of ALG 3
    double f_over_ell; double f_over_ell2;
    for (int ell = 0; ell < (2 * N) + 1; ell++) { //index starts from 0, so we stop when ell = 2N
        
        //firstly add f_over_ell
        if (i == N) { //we don't want to check X[i+1] in this case (index error), and in this case f_over_ell is always f(D') + x_j!
            check_if_i = false; //means the upcoming point f_over_ell is not a point of f(D) + x_i but rather f(D') + x_j
            f_over_ell = varphi + X[j + 1];
        }
        else if (j == N) { //same for X[j+1]
            check_if_i = true;
            f_over_ell = X[i + 1];
        }
        else { //otherwise check whether the next point is f(D) + x_i or f(D') + x_j
            if (X[i + 1] < varphi + X[j + 1]) { //means f(D) +  X[i+1] comes next,
                check_if_i = true;
                f_over_ell = X[i + 1];
            }
            else {
                check_if_i = false;
                f_over_ell = varphi + X[j + 1];
            }
        }
        length_ell = f_over_ell - f_under_ell; //Len(A_ell)

        //before the next loop starts, initialize some stuff!
        i2 = -1; // init
        j2 = -1; // init
        //same for varphi2 -- init
        if (varphi2 > 0) {//means f(D) < f(D') so f(D) + x_i comes first in the common refinement
            i2++;
            f_under_ell2 = X[i2];
        }
        else {
            j2++;
            f_under_ell2 = X[j2] + varphi2;
        }
        

        
        for (int ell2 = 0; ell2 < (2 * N) + 1; ell2++){
            if (i2 == N) { //we don't want to check X[i+1] in this case (index error), and in this case f_over_ell is always f(D') + x_j!
                check_if_i2 = false; //means the upcoming point f_over_ell is not a point of f(D) + x_i but rather f(D') + x_j
                f_over_ell2 = varphi2 + X[j2 + 1];
            }
            else if (j2 == N) { //same for X[j+1]
                check_if_i2 = true;
                f_over_ell2 = X[i2 + 1];
            }
            else { //otherwise check whether the next point is f(D) + x_i or f(D') + x_j
                if (X[i2 + 1] < varphi2 + X[j2 + 1]) { //means f(D) +  X[i+1] comes next,
                    check_if_i2 = true;
                    f_over_ell2 = X[i2 + 1];
                }
                else {
                    check_if_i2 = false;
                    f_over_ell2 = varphi2 + X[j2 + 1];
                }
            }
            //define length
            length_ell2 = f_over_ell2 - f_under_ell2; //Len(A_ell)
            //now check the three cases
            // case 1 (the first logical means that if the length of interval is zero don't add the probabilities since it will reduce sparsity of LP)
            if (abs(length_ell) <= pow(10, -6) || abs(length_ell2) <= pow(10, -6) || f_under_ell >= X[N] - pow(10, -6)  || f_over_ell <= X[0] + pow(10, -6) || f_under_ell2 >= X[N] - pow(10, -6)  || f_over_ell2 <= X[0] + pow(10, -6)) { //new addition at the end -- || i == p.size()-1 don't add this to the event if p_N comes in
                // do nothing: don't add this interval to the worst-case A \subseteq Omega
                indexD[ell][ell2] = -1; //means no access
                indexD2[ell][ell2] = -1; //means no access
                indexDD[ell][ell2] = -1; //means no access
                indexDD2[ell][ell2] = -1; //means no access
            }
            else if (f_under_ell >= X[N] + varphi - pow(10, -6) || f_over_ell <= X[0] + varphi + pow(10, -6) || f_under_ell2 >= X[N] + varphi2 - pow(10, -6) || f_over_ell2 <= X[0] + varphi2 + pow(10, -6)) { //means f(D') + X cannot access so always add!
                //add it!
                indexD[ell][ell2] = i; //means no access
                indexD2[ell][ell2] = i2; //means no access
                indexDD[ell][ell2] = -1; //means no access
                indexDD2[ell][ell2] = -1; //means no access
                violation = violation +  p[i][i2].get(GRB_DoubleAttr_X) ; //extra violation added by this A_\ell
            }
            else { //else check the condition
//                cout << i << i2 << j << j2 << "\n" << endl;
                condition = p[i][i2].get(GRB_DoubleAttr_X) - exp(epsilon)*(p[j][j2].get(GRB_DoubleAttr_X));
                if (condition > pow(10, -8)) {// add both in this case of violation
                    indexD[ell][ell2] = i; //means no access
                    indexD2[ell][ell2] = i2; //means no access
                    indexDD[ell][ell2] = j; //means no access
                    indexDD2[ell][ell2] = j2; //means no access
                    violation = violation + condition;
                }
                else { //don't add either
                    indexD[ell][ell2] = -1; //means no access
                    indexD2[ell][ell2] = -1; //means no access
                    indexDD[ell][ell2] = -1; //means no access
                    indexDD2[ell][ell2] = -1; //means no access
                }
            }
            
            f_under_ell2 = f_over_ell2; //the next step's A_ell will start from the end of this interval
            if (check_if_i2){
                i2++;
            }
            else{
                j2++;
            }
        }
        
        f_under_ell =f_over_ell; //the next step's A_ell will start from the end of this interval
        if (check_if_i) { //means that the next interval will start with an f(D) + x_i point so increase i
            i++;
        }
        else { //same for j
            j++;
        }
    }
}

// meta algorithm below
double optimization_alg(const vector<double>& X, double f_overline, vector<vector<double>>& p_decisions, int obj, int monoton) {//decisions is a dummy to keep optimized p vector
    auto start = high_resolution_clock::now(); //start counting the time
    vector<double> Delta = differences(X); remove(Delta); upper_lower(Delta, f_overline); sortVector(Delta);

    int N = (int)X.size() - 1;
    

    double violation;
    double worst_violation;              //will keep worst violation out of all violations in Delta


    
//    double index; double indexd; double index2; double indexd2;
    // ---- start GUROBI
    env = new GRBEnv(true); //start GUROBI environment
    //env->set("LogFile", "try_new_log.log"); //x-> operation is the same as *x.
    env->set(GRB_IntParam_OutputFlag, 0);   //mute GUROBI
    env -> set(GRB_IntParam_Method, 1); // DUAL SIMPLEX?

    //env->set(GRB_IntParam_Method, 1);
    env->start();                           //I start here, before the while loop, which will make the model warm-start all the time!
    model = new GRBModel(*env);             // Start
    
    for (auto& row : p) { //clear the rows
        row.clear(); // Clears the contents of each row, but keeps the empty rows
    }
    // Create variables
    for (int i = 0; i < N; ++i) {
        p.push_back(vector<GRBVar>(N));
        for (int j = 0; j < N; ++j) {
            p[i][j] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS); // Adjust bounds and variable type as necessary
        }
    }


    // Set objective: e.g., min amplitude (linear)!
    GRBLinExpr objchoice = 0;
    if (obj == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                objchoice += p[i][j] * ((abs(X[i + 1] + X[i]) + abs(X[j + 1] + X[j]))/2.0);
            }
        }
    }
    else if(obj == 2){
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                objchoice += p[i][j] * (min(abs(X[i]), abs(X[i + 1])) + min(abs(X[j]), abs(X[j + 1])));
            }
        }
    }
    model->setObjective(1 * objchoice, GRB_MINIMIZE);
    
    // Add constraint: sum(p) ==1 (most probably there is a shorter/faster implementation)
    GRBLinExpr sum_lhs = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum_lhs += p[i][j];
        }
    }
    model->addConstr(sum_lhs == 1.0); // both distributions sum to 1
    
    
    if (monoton == 1) { //add monotonicity constraints
        for (int j = 0; j < N; j++){
            for (int i = 0; i <= (N-2); i++) {
                if (X[i + 2] <= 0) {
                    model->addConstr(p[i][j] / (X[i + 1] - X[i]) - p[i + 1][j] / (X[i + 2] - X[i + 1]) <= 0);
                }
                else if (X[i] >= 0) {
                    model->addConstr(p[i][j] / (X[i + 1] - X[i]) - p[i + 1][j] / (X[i + 2] - X[i + 1]) >= 0);
                }
                else{
                    model->addConstr(p[i][j] / (X[i + 1] - X[i]) - p[i + 1][j] / (X[i + 2] - X[i + 1]) == 0);
                }
            }
        }
        for (int i = 0; i < N; i++){
            for (int j = 0; j <= (N-2); j++) {
                if (X[j + 2] <= 0) {
                    model->addConstr(p[i][j] / (X[j + 1] - X[j]) - p[i][j+1] / (X[j + 2] - X[j + 1]) <= 0);
                }
                else if (X[j] >= 0) {
                    model->addConstr(p[i][j] / (X[j + 1] - X[j]) - p[i][j+1] / (X[j + 2] - X[j + 1]) >= 0);
                }
                else{
                    model->addConstr(p[i][j] / (X[j + 1] - X[j]) - p[i][j+1] / (X[j + 2] - X[j + 1]) == 0);
                }
            }
        }
    }

    
    GRBLinExpr worst_case = 0;
    GRBLinExpr possible_worst_case = 0;
    int iterations = 0; //each iteration of the while loop
    while (true) { //iterate algorithm until no more violation
        iterations++;
        worst_violation = 0; //start with zero (this variable will keep the worst violation out of all the worst-case event of each varphi \in Delta
        worst_case = 0;
        
        model->update();
        model->optimize();
        // slack deletion below |
        //check all \varphi \in \Delta below
        for (int varphi = 0; varphi < Delta.size(); varphi++) { //check worst case of all varphi \in \Delta
            for (int varphi2 = 0; varphi2 < Delta.size(); varphi2++) { //check worst case of all varphi \in \Delta
                
                
                if (obj == 2){ //this needs to be obj = 2 to be fully correct
                    if ( pow(Delta[varphi],2) + pow(Delta[varphi2],2)  <= pow(f_overline,2) + pow(10,-8)  ){ // otherwise this does not fall in the ell_2 ball
//                        worst_event(epsilon, Delta[varphi], Delta[varphi2], X, p, indexD, indexD2,  indexDD,indexDD2, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
                        //                    cout << violation << endl;
                        if (violation >= worst_violation + pow(10, -2)) { //if violation found in worst_event (ALG3) is largest so far, update!
                            //                        cout << violation << endl;
                            //NOTE: Here I am checking every solution and taking the WORST violated one.
                            //save the worst-case scenario solution
//                            worst_D = indexD;
//                            worst_D2 = indexD2;
//                            worst_Dd = indexDD;
//                            worst_Dd2 = indexDD2;
                            worst_violation = violation;
                        }
                    }
                }
                else if (obj == 0){
                    bool checker_temp = false; //checker_temp will keep track whether previous stage was strictly included in the ball
                    double val1 = 100.0;
                    double val2 = 100.0;
                    if (Delta[varphi] > pow(10,-5)){
                        val1 = Delta[varphi - 1];
                    }
                    else if (Delta[varphi] < -pow(10,-5)){
                        val1 = Delta[varphi + 1];
                    }
                    else if (abs(Delta[varphi]) <= pow(10,-8)){
                        val1 = 0.0;
                    }
                    
                    if (Delta[varphi2] > pow(10,-5)){
                        val2 = Delta[varphi2 - 1];
                    }
                    else if (Delta[varphi2] < -pow(10,-5)){
                        val2 = Delta[varphi2 + 1];
                    }
                    else if (abs(Delta[varphi2]) <= pow(10,-8)){
                        val2 = 0.0;
                    }
                    
                    if( pow(val1,2) + pow(val2, 2) <= pow(f_overline,2) - pow(10,-8)){
                        checker_temp = true;
                    }
                        
                        
                    if ( pow(Delta[varphi],2) + pow(Delta[varphi2],2)  <= pow(f_overline,2) + pow(10, -8)  || checker_temp){ // otherwise this does not fall in the ell_2 ball
//                        worst_event(epsilon, Delta[varphi], Delta[varphi2], X, p, indexD, indexD2,  indexDD,indexDD2, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
                        worst_event_to_add(epsilon, Delta[varphi], Delta[varphi2], X,  p, violation, possible_worst_case);
                        if (violation >= worst_violation + pow(10, -7)) { //if violation found in worst_event (ALG3) is largest so far, update!
                            //                        cout << violation << endl;
                            //NOTE: Here I am checking every solution and taking the WORST violated one.
                            //save the worst-case scenario solution
                            worst_violation = violation;
                            worst_case = possible_worst_case;
                            
                        }
                    }
                    
                }
            }
        }
        //time to check whether the worst violation among all worst-case events is feasible or infeasible (stopping condition)
        if (worst_violation <= delta + pow(10,-8)) { //+0.00001 is for solver precision
            break; //no more violations
        }
        else { //otherwise add the worst constraint
//            add_dp_constr(sum_lhs, X, N, worst_D, worst_Dd,worst_D2, worst_Dd2, index, indexd, index2, indexd2);
            model->addConstr(worst_case <= delta);
        }
    }
    auto stop = high_resolution_clock::now(); //end time
    auto duration = duration_cast<milliseconds>(stop - start); //time passed
    //****things to return
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            p_decisions[i][j] = p[i][j].get(GRB_DoubleAttr_X);
        }
    }
    cout << "number of iterations: " << iterations << endl;
    cout << "time (ms): " << duration.count() << endl;
    return duration.count();
}


double noise_amplitude_lower(const vector<double>& X, const vector<vector<double>>& p) {
    double ampl = 0;
    int N = static_cast<int>(p.size()); // Assuming p is a square matrix for simplicity.
//    cout << p.size() << endl;
//    cout << p[1].size() << endl;
//    cout << X.size() << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ampl += p[i][j] * (min(abs(X[i]), abs(X[i + 1])) + min(abs(X[j]), abs(X[j + 1])));
        }
    }
    return ampl; // No division by 2 here as we're directly computing the desired sum
}

double noise_amplitude(const vector<double>& X, const vector<vector<double>>& p) {
    double ampl = 0;
    int N = static_cast<int>(p.size()); // Assuming p is a square matrix for simplicity.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ampl += p[i][j] * ((abs(X[i + 1] + X[i]) + abs(X[j + 1] + X[j])) / 2.0);
        }
    }
    return ampl; // No division by 2 here as we're directly computing the desired sum
}


double noise_power(const vector<double>& X, const vector<vector<double>>& p) {
    double power = 0;
    int N = static_cast<int>(X.size()) - 1; // Assuming X.size() gives N+1 points, and p is N x N.

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // For each cell, calculate the contribution to noise power.
            // Assuming a uniform distribution within each cell, we integrate x^2 over the area.
            double integral_x = (pow(X[i+1],3) - pow(X[i],3)) / 3.0;
            double integral_y = (pow(X[j+1],3) - pow(X[j],3)) / 3.0;
            power += p[i][j] * (pow(integral_x, 2) + pow(integral_y,2));
        }
    }

    return power;
}


/*
**************************************
Run and Test Functions
**************************************
*/
int main() { //this is the main 2D optimizing function....
    epsilon = 1.0; //desired epsilon
    delta =  0.1;  //desired delta
    
    double f_overline = 1.0; //sensitivity parameter in ell_2 norm!! (new)

    double increments= 1.0; //increments in each dimension
    
    //************************* Truncated Laplace
    double truncated_quantity = (exp(epsilon) - 1) / (2 * delta);
    double limits = 4; // this is not "that" simple. Here, because we worked a case where -+ 4 support is sufficient we fix that. For dynamically adjusting this, please use the practice as in the other codes.

    vector<double> X; arange(X, -limits, limits, increments); //input: construct the noise X -> here there is issue -- Also, here we are optimizing over a uniform partitioning noise vector.
//    print_vector(X);
    //************************* What to Optimize
    int objchoice = 0; //Objective Function: 0: amplitude, 1: power , 2: amplitude lb (you can add power lb, etc.)
    int monoton = 0; //0: nothing, 1: impose monotonicity -- latter is faster but does not give admissible lower bounds (that is, monotonic LB > general UB is possible)

    //************************* Truncated Stats
    //cout << "Gaussian Mechanism has amplitude = " << sqrt(2*log(1.25/delta)*pow(f_overline/epsilon, 2))*(sqrt(2)/sqrt(atan(1)*4)) << endl;
    vector<double> results;
    results.push_back((f_overline / epsilon) * log(1 + truncated_quantity));
    results.push_back((f_overline / epsilon) * (1 - (log(1 + truncated_quantity) / truncated_quantity)));
    results.push_back(2 * pow((f_overline / epsilon), 2) * (1 - ((log(1 + truncated_quantity) + 0.5 * pow(log(1 + truncated_quantity), 2)) / truncated_quantity)));
    cout << "Truncated Laplace has support [-A, A] where A= " << results[0] << endl;
    cout << "Truncated Laplace expected amplitude is: " << results[1] << endl;
    cout << "Truncated Laplace expected power is: " << results[2] << endl;
    //************************* Optimization Stats

    int N = (int)X.size() - 1;

    vector<vector<double>> p_decisions(N, vector<double>(N, 0.0)); //start N * N 2D probabilities
//    for (int i = 0; i < (int)p_decisions[0].size(); i++){
//        print_vector(p_decisions[i]);
//    }
//    double amount = 0.2;
    double duration = optimization_alg(X, f_overline, p_decisions, objchoice, monoton);
    results.push_back(noise_amplitude(X, p_decisions));
    results.push_back(noise_amplitude_lower(X, p_decisions));
//    results.push_back(noise_power(X, p_decisions));

//    print_vector(p_decisions[1]);

    cout << "Our amplitude is: " << results[3] << endl;

    cout << "Our amplitude lower bound is: " << results[4] << endl;
    
//    cout << "All for increments: " << results[5] << endl;
    results.push_back(duration/1000.0);
    cout << "Duration (s): " << results[6] << endl;
    //cout << "Here is the optimal solution:" << endl;
    //print_vector(p_decisions);

    return 0;
}



int main_worst() { //this is the worst-case finder. That is, if you give 2 distributions, we compute the product distribution's "delta" guarantees.
    epsilon = 5.0; //desired epsilon
    delta =  0.2;  //desired delta
    
    double f_overline =1.0; //sensitivity parameter in ell_2 norm!! (new)
    
    double max_eps = 5.0;
    vector<double> epsilon_grid = {};
    double step = round(max_eps / 20.0 * pow(10,6)) / pow(10,6);  // 19 intervals to get 20 elements
    epsilon_grid.clear();
    for (double ee = step; ee <= max_eps + pow(10, -6); ee += step){
        epsilon_grid.push_back(ee);
    }
    
    vector<double> delta_grid = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,0.08,0.09, 0.1,0.11,0.12,0.13,0.14, 0.15,0.16, 0.17, 0.18, 0.19, 0.2};
    vector<double> increments_grid = {0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0};
    
    double increments = increments_grid[0];
    
    double limits = 20.0; // because in our experiments we took limits = 20
    vector<double> X; arange(X, -limits, limits, increments); //input: construct the noise X -> here there is issue
    
    int N = (int)X.size() - 1;
    
    // Start these once only
    double worst_violation = 0.0;
    vector<double> Delta = differences(X); remove(Delta); upper_lower(Delta, f_overline); //set of potential worst-case [f(D') - f(D)] values, remove duplicates, remove the elements out of the range, add limits
    sortVector(Delta);
//    Delta = {Delta[0],0.0, Delta[Delta.size()-1]};
    
    //dummies to be used later
    vector<vector<double>> indexD(2 * N + 1, vector<double>(2 * N + 1));    //unassigned -- all intervals
    vector<vector<double>> indexD2(2 * N + 1, vector<double>(2 * N + 1));    //unassigned -- all intervals
    vector<vector<double>> indexDD(2 * N + 1, vector<double>(2 * N + 1));   //unassigned
    vector<vector<double>> indexDD2(2 * N + 1, vector<double>(2 * N + 1));   //unassigned


    //
    
    // FIRST DISTRIBUTION
    double best_amplitude = 100.0;
    
    for (int first_index = 0; first_index < 20; first_index++){
        for (int second_index = 0; second_index < 20; second_index++){
            if (epsilon_grid[first_index] + epsilon_grid[second_index] < epsilon - pow(10, -6) || delta_grid[first_index] + delta_grid[second_index] < delta - pow(10, -6)  ){
//                continue;
//                double xyz = 3.0;
                continue;
                // if yes, don't continue the next stuff, iterate the loop instead
            }
            double violation;
            epsilon = epsilon_grid[first_index]; //middle
            delta = delta_grid[first_index];
            std::stringstream stream; //start a stringstream
            stream << std::fixed << std::setprecision(3) << increments; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
            std::stringstream epsstream; //start a stringstream
            epsstream << std::fixed << std::setprecision(3) << epsilon;
            std::stringstream deltastream; //start a stringstream
            deltastream << std::fixed << std::setprecision(3) << delta;
            std::vector<double> results_1 = read_csv_to_vector("/Results/Rebuttal/inc-" + stream.str() + "-"+ epsstream.str()+ "-" + deltastream.str() + ".csv");
        //    print_vector(results_1);
            std::vector<double> dist_1 = read_csv_to_vector("/Results/Rebuttal/dist-inc-" + stream.str() + "-"+ epsstream.str()+ "-" + deltastream.str() + ".csv");
            // SECOND DISTRIBUTION
            epsilon = epsilon_grid[second_index]; //middle
            delta = delta_grid[second_index];
            std::stringstream epsstream2; //start a stringstream
            epsstream2 << std::fixed << std::setprecision(3) << epsilon;
            std::stringstream deltastream2; //start a stringstream
            deltastream2 << std::fixed << std::setprecision(3) << delta;
            std::vector<double> results_2 = read_csv_to_vector("/Results/Rebuttal/inc-" + stream.str() + "-"+ epsstream2.str()+ "-" + deltastream2.str() + ".csv");
        //    print_vector(results_2);
            std::vector<double> dist_2 = read_csv_to_vector("/Results/Rebuttal/dist-inc-" + stream.str() + "-"+ epsstream2.str()+ "-" + deltastream2.str() + ".csv");
            
            //construct the new p
            vector<vector<double>> p_decisions(N, vector<double>(N, 0.0)); //start N * N 2D probabilities
            for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    p_decisions[i][j] = dist_1[i]*dist_2[j];
                }
            }
            
            //amplitude compute
            double combined_amplitude = results_1[3] + results_2[3];
            
            //check if the worst violation <= 0.2
            epsilon = max_eps; // fix the epsilon first
            worst_violation = 0.0;
            bool cond_to_br = false;
            for (int varphi = 0; varphi < Delta.size(); varphi++) { //check worst case of all varphi \in \Delta
                for (int varphi2 = 0; varphi2 < varphi; varphi2++) { //check worst case of all varphi \in \Delta
                    if (cond_to_br == false){
                        bool checker_temp = false; //checker_temp will keep track whether previous stage was strictly included in the ball
                        double val1 = 100.0;
                        double val2 = 100.0;
                        if (Delta[varphi] > pow(10,-5)){
                            val1 = Delta[varphi - 1];
                        }
                        else if (Delta[varphi] < -pow(10,-5)){
                            val1 = Delta[varphi + 1];
                        }
                        else{
                            val1 = 0.0;
                        }
                        
                        if (Delta[varphi2] > pow(10,-5)){
                            val2 = Delta[varphi2 - 1];
                        }
                        else if (Delta[varphi2] < -pow(10,-5)){
                            val2 = Delta[varphi2 + 1];
                        }
                        else{
                            val2 = 0.0;
                        }
                        
                        if( pow(val1,2) + pow(val2, 2) <= pow(f_overline,2) - pow(10,-4)){
                            checker_temp = true;
                        }
                            
                    
                        if ( pow(Delta[varphi],2) + pow(Delta[varphi2],2)  <= pow(f_overline,2) + pow(10, -8) || checker_temp ){ // otherwise this does not fall in the ell_2 ball
//                            worst_event_static(epsilon, Delta[varphi], Delta[varphi2], X, p_decisions, indexD, indexD2,  indexDD,indexDD2, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
                            worst_event_improved(epsilon, Delta[varphi], Delta[varphi2], X, p_decisions, violation);
                            //                    cout << violation << endl;.
                            if (violation >= worst_violation + pow(10, -8)) { //if violation found in worst_event (ALG3) is largest so far, update!
                                //                        cout << violation << endl;
                                //NOTE: Here I am checking every solution and taking the WORST violated one.
                                //save the worst-case scenario solution
                                worst_violation = violation;
                                if (worst_violation > 0.2){
                                    cond_to_br = true;
                                }
                            }
                        }
                    }
                }
            }
            if (worst_violation <= 0.2 + pow(10,-3) && combined_amplitude <= best_amplitude - pow(10, -6) ){
                cout << "worst violation is:" << worst_violation << endl;
                cout << "corresponding amplitude is:" << combined_amplitude << endl;
                best_amplitude = combined_amplitude;
            }
        }
    }
    
    
    return 0;
}
