/*
Include Packages
*/
#include <iostream>     //input output
#include <vector>       //to add vector<type> variables
#include <chrono>       //to measure time ellapsed in a function
#include "gurobi_c++.h" //Gurobi's functions
#include <algorithm>    //I don't remember this -- delete?
#include <math.h>       //some maths functions
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;
using namespace std::chrono;
/*
Set Global Parameters (for example, same GUROBI model will be called again and again in the loops -- we need this)
*/
GRBEnv* env;            //GUROBI environment
GRBModel* model;        //model
int overline_k =  42;     //this will be overriden -- number of noise distributions
vector<vector<GRBVar>> p(overline_k);        //create a p vector, for probabilities, and the type of vector is gurobi var
//vector<vector<GRBVar>> z(overline_k);        //create a z vector, for probabilities, and the type of vector is gurobi var
GRBConstr* c;           //will be used to get constraints
double epsilon;         //problem constant  -- main
double delta;           //problem constant  -- main
//std::ofstream f_dist; //start the file stream
double sm = pow(10, -8);//small number
/*
PRELIMINARY Functions
*/
void write_to_csv(const vector<double>& vect, ofstream& f) {    //give vector and file, this will write a new line with comma separated values of thevector
    f << std::setprecision(10) << vect[0];
    for (int i = 1; i < vect.size(); i++) {
        f << ',' << std::setprecision(10) << vect[i];
    }
    f << '\n';
}
void print_vector(const vector<double>& v) {                    //prints a <double> vector
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << '\t';
    } cout << endl;
}
void print_vector_int(const vector<int>& v) {                    //prints an <int> vector
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << '\t';
    }
}
void arange(vector<double>& result, double start, double end, double increase) { //create a vector from 'start' to 'end' by 'increase' jumps
    result.clear();
    //result.reserve((end - start) / increase);
    for (double i = start; i <= end + sm; i = i + increase) {
        result.push_back(round(i * pow(10,9)) / pow(10,9));
    }
}
void remove(vector<double>& v){                                  //removes duplicates from a given vector -- O(nlog(n))
    sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
}
void upper_lower_single(vector<double>& v, const double threshold) { //removes elements > upper_f || < - upper_f, and adds these two limits
    v.erase(std::remove_if(
        v.begin(), v.end(),
        [&](const double& x) {
            return x >= threshold || x <= -threshold; // condition s.t. we should remove from vector
        }), v.end());
    //now add end_points
    v.push_back(threshold);
    v.push_back(-threshold);
}

void upper_lower(vector<double>& v, const double threshold_left, const double threshold_right) { //removes elements that are not in the given limits & adds the limit points (piecewise linear breakpoints)
    v.erase(std::remove_if(
        v.begin(), v.end(),
        [&](const double& x) { // this was const int...
            return (x > threshold_right - sm || x - sm < threshold_left);         // condition s.t. we should remove from vector -- pretty sure about - sm because we will add it anyway
        }), v.end());
    //now add end_points
    v.push_back(threshold_right);
    v.push_back(threshold_left);
    remove(v); // to remove duplocates
}
vector<double> differences(const vector<double>& v){ //{x_i - x_j \ : \ i,j = 1, ..., N+1}
    int size = static_cast<int>(v.size()); // int size = v.size(); //changed to prevent compiler issues
    vector<double> diff; diff.reserve(size * size); //difference vector to return
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            diff.push_back( round((v[i] - v[j])*pow(10,7))/pow(10,7) ); // you can just keep v_i - v_j if this gives error later on
        }
    }
    return diff;
}
void add_dp_constr(const vector<double>& X, int N, int k, int m, vector<int>& indexD, vector<int>& indexDD, vector<double>& intervals) { // same as the function for single-noise case, but has k and m for the first dimension of the noise
    int index, indexd;
    GRBLinExpr lhs = 0.0; //start building the lhs of constraint
    for (int i = 0; i < 2 * N + 1; i++) {//all intervals in the comref will be 2N+1
        index = indexD[i]; //index of i: [f(D) + x_i, f(D') + x_{i+1}] includes this interval
        if (index >= 0) {//else it means we cannot reach here, so don't add!
            lhs += intervals[i] * (p[k][index] / (X[index + 1] - X[index]));//summation over i part of the constraint -> Len(A_ell)*p_i / Len(f(D) + x_i, f(D) + x_{i+1})
        }
        indexd = indexDD[i]; //index of j: [f(D') + x_j, f(D') + x_{j+1}] includes this interval
        if (indexd >= 0) {
            lhs -= exp(epsilon) * intervals[i] * (p[m][indexd] / (X[indexd + 1] - X[indexd]));
        }
    }
    model->addConstr(lhs <= delta);
}
void domains(vector<double>& underline_I, vector<double>& overline_I, const double F_lower, const double F_upper, const int overline_k, const int method) {// divide X by overline_k intervals with method , F is function sensitivity
    if (method == 0) { // equally dividing -- more methods later
        double jump_size = (double) (F_upper - F_lower) / (1.0*overline_k);
        cout << "Jump size is: " << jump_size << endl;
        double temp = F_lower;
        for (int i = 0; i < overline_k; i++) {
            underline_I[i] = temp;
            temp += jump_size;
            overline_I[i] = temp;
        }
    }
    else { // this is experimental for now
        cout << "error." << endl;
    }
}
void domains_updated(vector<double>& underline_I, vector<double>& overline_I, const double F_lower, double F_upper, double F_truth, double increments, const int overline_k, const int method, const double f_overline) {// divide X by overline_k intervals with method , F is function sensitivity
    // some additional Phi partitioning methods
    if (method == 0) { // equally dividing -- more methods later
        domains_updated(underline_I, overline_I, F_lower, F_upper, F_truth, increments, overline_k, method, f_overline);
    }
    else if (method == 1){ // == 1 means +- increments around the middle one, and the rest is equal
        double jump_size_lower = (double) (F_truth - increments - F_lower) / ((overline_k - 1)/2);
        double temp_lower = F_lower;
        for (int i = 0; i < (int) (overline_k - 1)/2; i++) {
            underline_I[i] = temp_lower;
            temp_lower += jump_size_lower;
            overline_I[i] = temp_lower;
        }
        underline_I [(int) (overline_k - 1)/2] = F_truth - increments;
        overline_I [(int) (overline_k - 1)/2] = F_truth + increments;
        double jump_size_upper = (double) (F_upper - F_truth - increments) / ((overline_k - 1)/2);
        double temp_upper = F_truth + increments;
        for (int i = (int) (overline_k - 1)/2 + 1; i < overline_k; i++){
            underline_I[i] = temp_upper;
            temp_upper += jump_size_upper;
            overline_I[i] = temp_upper;
        }
//        cout << "Printing vector now " << endl;
//        print_vector(underline_I);
//        print_vector(overline_I);
    }
    else if (method == 2) { // == 2 means +- increments around the middle one, and one more the same, but the rest is the same
        double jump_size_lower = (double) (F_truth - 3*increments - F_lower) / ((overline_k - 3)/2);
        double temp_lower = F_lower;
        for (int i = 0; i < (int) (overline_k - 3)/2; i++) {
            underline_I[i] = temp_lower;
            temp_lower += jump_size_lower;
            overline_I[i] = temp_lower;
        }
        //now let's do the middle one
        for (int i = (int) (overline_k - 3)/2; i < (int) (overline_k + 1)/2 + 1; i++) {
            underline_I[i] = temp_lower;
            temp_lower += 2*increments; //length will be 2beta in the central 3 distributions
            overline_I[i] = temp_lower;
        }
        //finally last bit
        double jump_size_upper = (double) (F_upper - F_truth - 3*increments) / ((overline_k - 3)/2);
        double temp_upper = F_truth + 3*increments;
        for (int i = (int) (overline_k + 1)/2 + 1; i < overline_k; i++){
            underline_I[i] = temp_upper;
            temp_upper += jump_size_upper;
            overline_I[i] = temp_upper;
        }
//        cout << "Printing vector now " << endl;
//        print_vector(underline_I);
//        print_vector(overline_I);
    }
    else if (method == 3) { // == 3 means +- increments around the middle one, then two more intervals on both sides each with delta f / 2 so that the delta_f much neighbors are all protected
        double jump_size_lower = (double) (F_truth - increments - f_overline - F_lower) /  ((overline_k - 5.0)/2.0);
        double temp_lower = F_lower;
        for (int i = 0; i < (int) (overline_k - 5)/2; i++) {
            underline_I[i] = temp_lower;
            temp_lower += jump_size_lower;
            overline_I[i] = temp_lower;
        }
        underline_I [(int) (overline_k - 5)/2] = F_truth - increments - f_overline;
        overline_I [(int) (overline_k - 5)/2] = F_truth - increments - (f_overline/2.0);
        underline_I [(int) (overline_k - 3)/2] = F_truth - increments - (f_overline/2.0);
        overline_I [(int) (overline_k - 3)/2] = F_truth - increments;
        underline_I [(int) (overline_k - 1)/2] = F_truth - increments;
        overline_I [(int) (overline_k - 1)/2] = F_truth + increments;
        underline_I [(int) (overline_k + 1)/2] = F_truth + increments;
        overline_I [(int) (overline_k + 1)/2] = F_truth + increments + (f_overline/2.0);
        underline_I [(int) (overline_k + 3)/2] = F_truth + increments + (f_overline/2.0);
        overline_I [(int) (overline_k + 3)/2] = F_truth + increments + (f_overline);
        //finally last bit
        double jump_size_upper = (double) (F_upper - F_truth - increments - f_overline) / ((overline_k - 5.0)/2.0);
        double temp_upper = F_truth + increments + (f_overline);
        for (int i = (int) (overline_k + 3)/2 + 1; i < overline_k; i++){
            underline_I[i] = temp_upper;
            temp_upper += jump_size_upper;
            overline_I[i] = temp_upper;
        }
        cout << "Printing vector now " << endl;
        print_vector(underline_I);
        print_vector(overline_I);
    }
}
void ranges_nonunif(vector<vector<vector<double>>>& Deltas, vector<vector<int>>& validity, const double f_overline, vector<double>& underline_I, vector<double>& overline_I, const vector<double>& X) {
    vector<double> Delta; //will be re-used
    vector<double> Diffs = differences(X); remove(Diffs); //prepare the differences x_i - x_j here as all Deltas will start from this
    double lhs, rhs;
    for (int k = 0; k < overline_k; k++) {
        for (int m = 0; m < overline_k; m++) {
            lhs = underline_I[m] - overline_I[k];
            rhs = overline_I[m] - underline_I[k];
        
            Delta = Diffs;
            Delta.erase(std::remove_if(Delta.begin(), Delta.end(), [&](const double& x) {
                return (x > f_overline + sm || x + sm < -f_overline);
                }), Delta.end());
            Delta.erase(std::remove_if(Delta.begin(), Delta.end(), [&](const double& x) {
                return (x > rhs - sm || x - sm < lhs);
                }), Delta.end());
            Deltas[k][m] = Delta;
            if (Delta.size() < 1){
                validity[k][m] = -1;
            }
            else{
                validity[k][m] = 1;
            }
//            cout << "k to m " << to_string(k) << " to " << to_string(m) << ": " <<endl;
//            cout << "lhs is " << lhs << endl;
//            cout << "rhs is " << rhs << endl;
//            cout << "f_overline is " << f_overline << endl;
//            print_vector(Delta);
//            cout << endl;
        }
    }
}
void ranges_old(vector<vector<vector<double>>>& Deltas, vector<vector<int>>& validity, const double f_overline, vector<double>& underline_I, vector<double>& overline_I, const vector<double>& X) {
    vector<double> Delta; //will be re-used
    cout << f_overline << endl;
    vector<double> Diffs = differences(X); remove(Diffs); upper_lower_single(Diffs, f_overline); //prepare the differences x_i - x_j here as all Deltas will start from this
    double lhs, rhs;
    for (int k = 0; k < overline_k; k++) {
        for (int m = 0; m < overline_k; m++) {
            lhs = round((underline_I[m] - overline_I[k])*10000.0)/10000.0;
            rhs = round((overline_I[m] - underline_I[k])*10000.0)/10000.0;
            Delta = Diffs;
//
//
//            Delta.erase(std::remove_if(Delta.begin(), Delta.end(), [&](const double& x) {
//                return (x > rhs - sm || x - sm < lhs);
//                }), Delta.end());
            
            bool deleted_rhs = false;
            bool deleted_lhs = false;
            
            Delta.erase(std::remove_if(Delta.begin(), Delta.end(), [&](const double& x) {
                if (x > rhs - sm) {
                    if (rhs <= f_overline + sm && rhs + sm >= - f_overline){
                        deleted_rhs = true;
                    }
                    return true;
                }
                else if (x < lhs + sm){
                    if (lhs + sm >= -f_overline && lhs <= f_overline + sm){
                        deleted_lhs = true;
                    }
                    return true;
                }
                else {
                    return false;
                }
            }), Delta.end());
            
            if (deleted_rhs == true){
                Delta.push_back(rhs);
            }
            if (deleted_lhs == true){
                Delta.push_back(lhs);
            }
            
            Deltas[k][m] = Delta;
            if (Delta.size() < 1){
                validity[k][m] = -1;
            }
            else{
                validity[k][m] = 1;
            }
//                        cout << "k to m " << to_string(k) << " to " << to_string(m) << ": " <<endl;
//                        cout << "lhs is " << lhs << endl;
//                        cout << "rhs is " << rhs << endl;
//                        cout << "f_overline is " << f_overline << endl;
//                        std::sort(Delta.begin(), Delta.end());
//                        print_vector(Delta);
//                        cout << endl;
        }
    }
}

void ranges(vector<vector<vector<double>>>& Deltas, vector<vector<int>>& validity, const double f_overline, vector<double>& underline_I, vector<double>& overline_I, const vector<double>& X) { // for UB
    vector<double> Delta; //will be re-used
    vector<double> Diffs = differences(X); remove(Diffs); //prepare the differences x_i - x_j here as all Deltas will start from this
//    cout << "There u goooooo hon:" << endl;
//    print_vector(Diffs);
    double lhs, rhs;
    for (int k = 0; k < overline_k; k++) {
        for (int m = 0; m < overline_k; m++) {
            lhs = max(round((underline_I[m] - overline_I[k])*100000.0)/100000.0 , -f_overline);
            rhs = min(round((overline_I[m] - underline_I[k])*100000.0)/100000.0 , f_overline);
            
            if (rhs - lhs <= pow(10,-6)){ //means the intersection is a singleton -- discard
                validity[k][m] = -1; // no access
                Delta = {};
                Deltas[k][m] = Delta;
            }
            else{ //there is a non-zero length intersection -- time to iterate in it and collect varphi values
                Delta = Diffs; // copy the "Diffs" which is the {(pi_j - pi_j') * beta} collection
                upper_lower(Delta, lhs, rhs); //remove <= lhs and >= rhs and then add these
                Deltas[k][m] = Delta;
                
                if (Delta.size() < 1){
                    validity[k][m] = -1;
                }
                else{
                    validity[k][m] = 1;
                }
            }
            cout << "k to m " << to_string(k) << " to " << to_string(m) << ": " <<endl;
            cout << "lhs is " << lhs << endl;
            cout << "rhs is " << rhs << endl;
            cout << "f_overline is " << f_overline << endl;
            std::sort(Delta.begin(), Delta.end());
            print_vector(Delta);
            cout << endl;
        }
    }
}

void ranges_LB(vector<vector<vector<double>>>& Deltas, vector<vector<int>>& validity, const double f_overline, vector<double>& underline_I, vector<double>& overline_I, const vector<double>& X) { //use this version when |Phi_k(\beta)| and |I_i(\beta)| are identical. Otherwise use the other "ranges" one. This is the LB breakpoints for varphi
    vector<double> Delta; //will be re-used
    vector<double> Diffs = differences(X); remove(Diffs); //prepare the differences x_i - x_j here as all Deltas will start from this
    for (int k = 0; k < overline_k; k++) {
        for (int m = 0; m < overline_k; m++) {
            double quant = (m-k)*(X[1] - X[0]);
            Deltas[k][m] = { quant };
            if (abs((m-k)*(X[1] - X[0])) > f_overline + sm){
                validity[k][m] = -1;
            }
            else{
                cout << "k to m " << to_string(k) << " to " << to_string(m) << ": " <<endl;
                cout << quant << endl;
                cout << endl;
                validity[k][m] = 1;
            }
        }
    }
}

/*
ESSENTIAL Functions
*/
void greedy_worst_event(const double epsilon, const double varphi, const int k, const int m, vector<double>& X, const vector<vector<GRBVar>>& p, vector<int>& indexD, vector<int>& indexDD, vector<double>& intervals, double& violation) { //returns the worst-case violation for neighbors k to m for a fixed varphi
    
    //assume wlog that f(D) = 0 and f(D') = varphi
    double condition; //condition to check
    double length_ell; //length of each interval
    bool check_if_i; //checks if the interval starts with an f(D)+i point or not
    violation = 0; //initially no violation
    // // int N = X.size() - 1; // -1 because X has size N+1 and probability vector has size N (N intervals)
    int N = static_cast<int>(X.size()) - 1;
    int i = -1, j = -1; //position of the last noise index f(D) + x_i or f(D') + x_j. =-1 means no access or don't include this interval
    double f_under_ell; //this will keep the beginning of A_ell = [f_under_ell, f_over_ell]
    //****** Initialize the first f_underline_ell
    if (varphi - sm >= 0) {//means f(D) < f(D') so f(D) + x_i comes first in the common refinement
        i++;
        f_under_ell = X[i];
    }
    else {
        j++;
        f_under_ell = X[j] + varphi;
    }
    //****** main for loop of ALG 3
    double f_over_ell;
    for (int ell = 0; ell < (2 * N) + 1; ell++) { //index starts from 0, so we stop when ell = 2N
        //firstly add f_over_ell
        if (i == N) { //we don't want to check X[i+1] in this case  (index error), and in this case f_over_ell is always f(D') + x_j!
            check_if_i = false; //means the upcoming point f_over_ell is not a point of f(D) + x_i but rather f(D') + x_j
            f_over_ell = varphi + X[j + 1];
        }
        else if (j == N) { //same for X[j+1]
            check_if_i = true;
            f_over_ell = X[i + 1];
        }
        else { //otherwise check whether the next point is f(D) + x_i or f(D') + x_j
//            cout << ell << endl;
            if (X[i + 1] < varphi + X[j + 1] - sm) { //means f(D) +  X[i+1] comes next,
                check_if_i = true;
                f_over_ell = X[i + 1];
            }
            else {
                check_if_i = false;
                f_over_ell = varphi + X[j + 1];
            }
        }
        //define length
        length_ell = f_over_ell - f_under_ell; //Len(A_ell)
        intervals[ell] = length_ell; //We will keep the interval lengths since these will be used in the privacy constraints (length * probability)
        //now check the three cases
            // case 1 (the first logical means that if the length of interval is zero don't add the probabilities since it will reduce sparsity of LP)
        if (abs(length_ell) <= sm || f_under_ell + sm >= X[N] || f_over_ell <= X[0] + sm) { //after the first 3, the new conditions are for LB ->  || i <= 2 || i >= N-2
            // do nothing: don't add this interval to the worst-case A \subseteq Omega
            indexD[ell] = -1;   //-1 means no access
            indexDD[ell] = -1;  //-1 means no access
        }
        else if (f_under_ell + sm >= X[N] + varphi || f_over_ell <= X[0] + varphi + sm) { //means f(D') + X cannot access so always add!
            //add it!
            indexD[ell] = i;
            indexDD[ell] = -1;//-1 since it cannot access
            violation = violation + (length_ell * (p[k][i].get(GRB_DoubleAttr_X) / (X[i + 1] - X[i]))); //extra violation added by this A_\ell
        }
        else { //else check the condition
            condition = (p[k][i].get(GRB_DoubleAttr_X) / (X[i + 1] - X[i])) - (exp(epsilon) * (p[m][j].get(GRB_DoubleAttr_X) / (X[j + 1] - X[j]))); //condition to check
            if (condition >= sm) { //then STRICT violation, so add this interval!
                indexD[ell] = i;
                indexDD[ell] = j;
                violation = violation + (length_ell * condition);
            }
            else {                  //no violation so do not add this
                indexD[ell] = -1;
                indexDD[ell] = -1;
            }
        }
        f_under_ell = f_over_ell; //the next step's A_ell will start from the end of this interval
        if (check_if_i) { //means that the next interval will start with an f(D) + x_i point so increase i
            i++;
        }
        else { //same for j
            j++;
        }
    }
}
double optimization_alg(vector<double>& X, double F_lower, double F_upper, double f_overline, double increments, int overline_k, vector<vector<double>>& p_decisions, int option, double threshold, int obj, int break_method, vector<double>& obj_weight, int monoton) {
    //decisions is a dummy to keep optimized p vector
    auto start = high_resolution_clock::now(); //start counting the time
    vector<double> underline_I(overline_k), overline_I(overline_k); //these will be filled soon]
    if(break_method == 0){ // old one where we divide equally
        double jumps = round(((F_upper - F_lower)/overline_k)*1000000.0) / 1000000.0;
        arange(underline_I, F_lower, F_upper - pow(10,-6), jumps);
        arange(overline_I, F_lower + jumps, F_upper + pow(10, -6), jumps);
        print_vector(underline_I);
        print_vector(overline_I);
    }
    else if(break_method == 7){ // a manual break method example for the Phi partitioning
//        domains(underline_I, overline_I, F_lower, F_upper, overline_k, break_method); //fill underline, overline I for each overline_k region where we add different noises
        arange(underline_I, -20.0, 20.51, 1.0);
        double myconstant = 0.5;
        std::transform(underline_I.begin(), underline_I.end(), underline_I.begin(), [&myconstant](auto& c){return c*myconstant;});
        overline_I = std::vector<double>(underline_I.begin() + 1, underline_I.end());
//        overline_I.push_back(10.0);
//  limits
//        underline_I.push_back(10.0);
        overline_I.push_back(12.0);
        underline_I.insert(underline_I.begin(), -12.0);
        overline_I.insert(overline_I.begin(), -10.0);
        print_vector(underline_I);
        print_vector(overline_I);
        for (int i = 0; i < underline_I.size(); i++) {
            underline_I[i] = round(underline_I[i] * 1000000.0) / 1000000.0;
            overline_I[i] = round(overline_I[i] * 1000000.0) / 1000000.0;
          }
    }
    else{
        domains_updated(underline_I, overline_I, F_lower, F_upper, (F_lower + F_upper)/2, increments, overline_k, break_method, f_overline); //fill underline, overline I for each overline_k region where we add different noises
    }
    vector<vector<vector<double>>> Deltas(overline_k, vector<vector<double>>(overline_k, vector<double>(1))); //Deltas will be deltas for neigbors k-m where varphi\in Deltas in worst-case
    vector<vector<int>> validity(overline_k, vector<int>(overline_k)); //validity = -1 means no need to consider that k-m combination
    ranges(Deltas, validity, f_overline, underline_I, overline_I, X);
    //int N = X.size() - 1;
    int N = static_cast<int>(X.size()) - 1;
    vector<double> Delta;
    //dummies to be used later
    vector<int> indexD((2 * N) + 1);    //unassigned
    vector<int> indexDD((2 * N) + 1);   //unassigned
    vector<double> intervals((2 * N) + 1); //unassigned
    double violation; double worst_violation;              //will keep worst violation out of all violations in Delta
    vector<int> worst_D((2 * N) + 1);
    vector<int> worst_Dd((2 * N) + 1);
    vector<double> worst_intervals((2 * N) + 1);
    double slack;                       //will keep the slack for each constraints
    // ---- start GUROBI
    env = new GRBEnv(true); //start GUROBI environment
//    env->set("LogFile", "/Users/.../Desktop/try.log"); //x-> operation is the same as *x.
    env->set(GRB_IntParam_OutputFlag, 0);   //mute GUROBI
//    env->set(GRB_IntParam_Method, 1);
    env->start();                           //I start here, before the while loop, which will make the model warm-start all the time!
    model = new GRBModel(*env);             // Start
//    p.clear();                              //clear solution from the other iterations
    // Create variables
    //p = model.
    for (int k = 0; k < overline_k; k++) {
        for (int i = 0; i < N; i++) {//create p[0]... p[n-1]  (in total n)
            p[k].push_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS)); //p is between 0-1
        }
    }
    // Set objective: e.g., min amplitude (linear)!
    GRBLinExpr objchoice = 0.0;
    if (obj == 0) {
        for (int k = 0; k < overline_k; k++) {
            for (int i = 0; i < N; i++) {
                objchoice += obj_weight[k] * p[k][i] * abs( (X[i + 1] + X[i]))/2.0;
            }
        } //objchoice = objchoice/overline_k;
    } //E[|X|] where we assume all regions "k" have the same likelihood
    else if (obj == 1) { //power
        for (int k = 0; k < overline_k; k++) {
            for (int i = 0; i < N; i++) {
                objchoice += obj_weight[k] * p[k][i] * (pow(X[i + 1], 3) - pow(X[i], 3)) / (3 * (X[i + 1] - X[i]));
            }
        }
        //objchoice = objchoice / overline_k;
    }
    else if (obj == 2) { //lower bound optimization
        for (int k = 0; k < overline_k; k++) {
            for (int i = 0; i < N; i++) {
                objchoice += obj_weight[k] * (p[k][i] * min(abs(X[i]), abs(X[i + 1])));
            }
        }
    }
    
    else if (obj == 10) { //to the power 10
        for (int k = 0; k < overline_k; k++) {
            for (int i = 0; i < N; i++) {
                objchoice += obj_weight[k] * pow(X[i + 1], 10);
            }
        }
    }
    // charge opposite sign more -- prevents sign flips -- used in the PCD example.
    else if (obj == 11) { //based on the sign
        for (int k = 0; k < overline_k; k++) {
            if (overline_I[k] <= 0.0001){ // means that this interval is on the left of zero, so we will charge the positive side more
                for (int i = 0; i < N; i++) {
                    if(X[i+1] <= 0.0001){ //it can go this much
                        objchoice += obj_weight[k] * p[k][i] * abs( (X[i + 1] + X[i]))/2.0;
                    }
                    else{
                        objchoice += obj_weight[k] * p[k][i] * 10000 * abs( (X[i + 1] + X[i]))/2.0; // 10 times more cost for the positive side
                    }
                    
                }
            }
            else{
                for (int i = 0; i < N; i++) {
                    if(X[i] <= 0.0001){ //it can go this much
                        objchoice += obj_weight[k] * p[k][i] * 10000 * abs( (X[i + 1] + X[i]))/2.0;
                    }
                    else{
                        objchoice += obj_weight[k] * p[k][i] * abs( (X[i + 1] + X[i]))/2.0; // 10 times more cost for the positive side
                    }
                    
                }
            }
            
        }
    }
    // charge opposite sign more when it flips the sign
    else if (obj == 12) { //based on the sign
        for (int k = 0; k < overline_k; k++) {
            if (overline_I[k] <= 0.0001){ // means that this interval is on the left of zero, so we will charge the positive side more
                for (int i = 0; i < N; i++) {
                    if(underline_I[k] + X[i+1] <= 1.00){ //this means that in the worst case of X, all of the interval will still fall on the good side
                        objchoice += obj_weight[k] * p[k][i] * abs( (X[i + 1] + X[i]))/2.0;
                    }
                    else{
                        objchoice += obj_weight[k] * p[k][i] * 1000 * abs( (X[i + 1] + X[i]))/2.0; // 10 times more cost for the positive side
                    }
                    
                }
            }
            else{
                for (int i = 0; i < N; i++) {
                    if(overline_I[k] + X[i] <= -1.00){ //it can go this much
                        objchoice += obj_weight[k] * p[k][i] * 1000 * abs( (X[i + 1] + X[i]))/2.0;
                    }
                    else{
                        objchoice += obj_weight[k] * p[k][i] * abs( (X[i + 1] + X[i]))/2.0; // 10 times more cost for the positive side
                    }
                    
                }
            }
            
        }
    }
    
    model->setObjective(1 * objchoice, GRB_MINIMIZE);
    // Add constraint: sum(p) ==1 for all \overline{k} distributions (most probably there is a quicker implementation)
    GRBLinExpr sum_lhs;
    for (int k = 0; k < overline_k; k++) {
        sum_lhs = 0.0;
        for (int i = 0; i < N; i++) {
            sum_lhs += p[k][i];
        }
        model->addConstr(sum_lhs == 1.0);
    }
    
    //new UB constraints //delete if you don't want to control the worst-case loss
//    for (int k = 0; k < overline_k; k++) {
//        GRBLinExpr ub_each;
//        for (int i = 0; i < N; i++) {
//            ub_each += p[k][i] * abs( (X[i + 1] + X[i]))/2.0;
//        }
//        model->addConstr(ub_each <= 1.42);
//    }
    
    // force all the distributions to be same -- uncomment if  you want them (but this will give the data independent noise)
//    for (int k = 0; k < overline_k - 1; k++) {
//        for (int i = 0; i < N; i++) {
//            model->addConstr((p[k][i] - p[k+1][i]) == 0);
//        }
//    }
    // force all distributions to be symmetric
//    for (int k = 0; k < overline_k; k++) {
//        for (int i = 0; i < floor(N/2); i++) {
//            model->addConstr((p[k][i] - p[k][N-i-1]) == 0);
//        }
//    }
//
//    LHS of zero and RHS of zero cannot differ that much for any probability distribution
//    for (int k = 0; k < overline_k; k++) {
//        GRBLinExpr left;
//        GRBLinExpr right;
//        for (int i = 0; i < int (N/2); i++) {
//            left += p[k][i];
//        }
//        for (int i = int (N/2); i < N; i++) {
//            right += p[k][i];
//        }
//        model->addConstr(left - right <= 0.35);
//        model->addConstr(right - left <= 0.35);
//    }

    
    
    //monotonicity constraints //activate below if you want cross distributional monotonicity constraints
//    if (monoton == 1) {
//        for (int k = 0; k < overline_k; k++) {
//            for (int i = 0; i < N; i++) {//create p[0]... p[n-1]  (in total n)
//                z[k].push_back(model->addVar(0.0, 1.0, 0.0, GRB_BINARY)); //z is binary
//            }
//        }
//        for (int k = 0; k < overline_k; k++) {
//            for (int i = 0; i < N - 1; i++) {
//                model->addConstr(z[k][i] - z[k][i + 1] <= 0);
//                model->addConstr(p[k][i + 1] - p[k][i] + z[k][i] >= 0);
//                model->addConstr(p[k][i + 1] - p[k][i] - (1 - z[k][i]) <= 0);
//            }
//        }
//    }
    //monotonicity constraints-- simple intra-distributional monotonicity
    if (monoton == 1) {
        for (int k = 0; k < overline_k; k++){
            for (int i = 0; i <= (N-2); i++) {
                if (X[i + 2] <= 0) {
                    model->addConstr((p[k][i] / (X[i + 1] - X[i])) - (p[k][i + 1] / (X[i + 2] - X[i + 1])) <= 0.0);
                }
                else if (X[i] >= 0) {
                    model->addConstr((p[k][i] / (X[i + 1] - X[i])) - (p[k][i + 1] / (X[i + 2] - X[i + 1])) >= 0.0);
                }
                else{
                    model->addConstr( (p[k][i] / (X[i + 1] - X[i])) - (p[k][i + 1] / (X[i + 2] - X[i + 1])) == 0.0);
                }
            }
        }
    }
    else if (monoton == 2) { // impose monotonicity except for the middle
        for (int k = 0; k < overline_k; k++){
            if (k != (int) (overline_k - 1)/2){
                for (int i = 0; i <= (N-2); i++) {
                    if (X[i + 2] <= 0) {
                        model->addConstr(p[k][i] / (X[i + 1] - X[i]) - p[k][i + 1] / (X[i + 2] - X[i + 1]) <= 0);
                    }
                    else if (X[i] >= 0) {
                        model->addConstr(p[k][i] / (X[i + 1] - X[i]) - p[k][i + 1] / (X[i + 2] - X[i + 1]) >= 0);
                    }
                }
            }
        }
    }
    //rest of the code
    int iterations = 0; //each iteration of the while loop
    int stopping_condition = 1;
    while (stopping_condition == 1) { //iterate algorithm until no more violation
        iterations++;
        model->update();
        model->optimize();
        if (model->get(GRB_IntAttr_Status) != 2) {
            cout << "ERROR!!!!! Status: "<< model->get(GRB_IntAttr_Status) << endl;
        }
        // slack deletion below |
        if (iterations % option == 0) { //the if condition ensures that slack deletion will be implemented every 'option' steps
            cout << "deleting" << endl;
            c = model->getConstrs();
            for (int constiter = 0; constiter < model->get(GRB_IntAttr_NumConstrs); ++constiter) {
                slack = c[constiter].get(GRB_DoubleAttr_Slack);
                if (slack > threshold) {
                    model->remove(c[constiter]);
                }
            }
        }
        stopping_condition = 0; //unless we add constraints below in the for loop-- 0 means stop.
        int worst_k = -1, worst_m =-1; //just to prevent warnings
        for (int k = 0; k < overline_k; k++) {
            for (int m = 0; m < overline_k; m++) {
                if (validity[k][m] == 1) { //otherwise we set these are not valid (cannot access)
                    worst_violation = 0.0; //start with zero (this variable will keep the worst violation out of all the worst-case event of each varphi \in Delta
                    //check all \varphi \in \Delta below
                    Delta = Deltas[k][m]; //deltas to look for
//                    cout << "k to m: " << to_string(k) << " and " << to_string(m) << endl;
//                    print_vector(Delta);
                    for (int varphi = 0; varphi < Delta.size(); varphi++) { //check worst case of all varphi \in \Delta
                        greedy_worst_event(epsilon, Delta[varphi], k, m, X, p, indexD, indexDD, intervals, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
//                        greedy_worst_event(epsilon, round(Delta[varphi] * 1000.0) / 1000.0, k, m, X, p, indexD, indexDD, intervals, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
                        if (violation > worst_violation + pow(10,-8)) {
                            worst_D = indexD;
                            worst_Dd = indexDD;
                            worst_intervals = intervals;
                            worst_violation = violation;
                            worst_k = k;
                            worst_m = m;
                        }
                    }
                    //time to check whether the worst violation among all worst-case events is feasible or infeasible (stopping condition)
                    if (worst_violation > delta + pow(10,-8)) { //add the worst constraint if it is infeasible
                        //cout << worst_violation << endl;
                        stopping_condition = 1; //we added at least one constraint so we should keep optimizing
//                        cout << "worst viol: " << worst_violation << endl;
//                        cout << "optimal objective value so far: " << objchoice.getValue() << endl;
                        add_dp_constr(X, N, worst_k, worst_m, worst_D, worst_Dd, worst_intervals);
                    }
                }
            }
        }
    }
    auto stop = high_resolution_clock::now(); //end time
    auto duration = duration_cast<milliseconds>(stop - start); //time passed
    //****things to return
    for (int k = 0; k < overline_k; k++) {
        for (int i = 0; i < N; i++) {//save all the GUROBI soltuions of p to the vector 'decisions'
            p_decisions[k][i] = p[k][i].get(GRB_DoubleAttr_X);
        }
    }
    cout << "number of iterations: " << iterations << endl;
    cout << "time (ms): " << duration.count() << endl;
    cout << "optimal objective value: " << objchoice.getValue() << endl;
    return duration.count();
}
double noise_amplitude(vector<double> X,vector<double> obj_weight, vector<vector<double>> p_decisions) {
    double ampl = 0.0;
    for (int k = 0; k < overline_k; k++) {
        for (int i = 0; i < p_decisions[0].size(); i++) {
            ampl += obj_weight[k] * p_decisions[k][i] * abs((X[i] + X[i + 1]))/2;
        }
    }
    return ampl;
}
double noise_amplitude_lower(vector<double> X,vector<double> obj_weight, vector<vector<double>> p_decisions) {
    double ampl = 0.0;
    double expect = 0.0;
    for (int k = 0; k < overline_k; k++) {
        expect = 0.0;
        for (int i = 0; i < p_decisions[0].size(); i++) {
            //cout << obj_weight[k] << endl;
            expect += (p_decisions[k][i] * min(abs(X[i]), abs(X[i + 1])));
        }
        ampl += obj_weight[k]*expect;
    }
    return ampl;
}
/*
Main function
*/
int main_multiple(double f_overline, double F_lower, double F_upper, int tr, double increments, int job_nr_to_read) {
    epsilon = 1.0; //desired epsilon
    delta = 0.2; //desired delta
    int break_method = 0; //split to four regions, break regions equally.
    //************************* Truncated Laplace
    double truncated_quantity = (exp(epsilon) - 1) / (2 * delta);
    //double limit_of_trunc = (f_overline / epsilon) * log(1 + truncated_quantity);
    // figure out parameters
    //************************* Optimizing DP
//    double limits = ceil(limit_of_trunc); //noise varies from -limits to limits
    cout << "there will be: " << overline_k << " distributions to optimize." <<endl;
//    cout << "limit taken initially: " << limits << endl;
//    double increments = 2*limits/number_of_variables; //aka beta -- take it so that overall we have "number_of_variables" decision variables
//    cout << "beta taken: " << increments << endl;
//    increments = f_overline/floor(f_overline/increments); //round the "increments" so that "f_overline" is divisible (we also state in the paper that Delta f / beta is divisible)
//    cout << "beta fixed for f_overline divisibility: " << increments << endl;
//    limits = number_of_variables*increments / 2; //same rounding for [-Lim, +Lim] (in papers notation [-Lbeta, Lbeta]
//    cout << "Corrected limits: " << limits << endl;
    //below just to update upper neighboring limits
    double F_lower_updated = F_lower;
    double F_upper_updated = F_upper;
    //below are some example "corrections" for the cases for the corection that needs to be taken for the assumptions to be satisfied
//    double current_length_upper, desired_upper_frac, desired_upper, F_upper_updated, current_length_lower, desired_lower_frac, desired_lower, F_lower_updated;
//    if (break_method == 0){
//        cout << "Now, let us fix function limits to be divisible by beta. Originally F_lower: " << F_lower << " and F_upper: " << F_upper << endl;
//        double desired_length = ceil((F_upper - F_lower) / (1.0 *overline_k));
//        F_upper_updated = desired_length*overline_k + F_lower;
//        F_lower_updated = F_lower;
//        cout << "UNIFORM SPLIT -> F_upper updated as: " << F_upper_updated  <<endl;
//    }
//    else if (break_method == 1){
//        cout << "Now, let us fix function limits to be divisible by beta. Originally F_lower: " << F_lower << " and F_upper: " << F_upper << endl;
//        current_length_upper = (F_upper - F_truth - increments)/(0.5*(overline_k - 1)); //length of each piece in the upper half from the true value
//        desired_upper_frac = ceil(current_length_upper / increments); //desired length fraction of each piece
//        desired_upper = desired_upper_frac * increments; //desired length of each upper piece
//        cout << "The lx§ength of one upper interval is thus " << current_length_upper << " but desired length is: " << desired_upper << endl;
//        F_upper_updated = ((double) desired_upper*(0.5 * (overline_k - 1)) + (F_truth + increments));
//        cout << "Hence, F_upper is updated as " << F_upper_updated << endl;
//        //below just to update lower neighboring limits
//        current_length_lower = (F_truth - increments - F_lower)/(0.5*(overline_k - 1)); //length of each piece in the upper half from the true value
//        desired_lower_frac = ceil(current_length_lower / increments); //desired length fraction of each piece
//        desired_lower = desired_lower_frac * increments; //desired length of each upper piece
//        cout << "The length of one lower interval is thus " << current_length_lower << " but desired length is: " << desired_lower << endl;
//        F_lower_updated = ((double) (F_truth - increments) - desired_lower*(0.5 * (overline_k - 1)) );
//        cout << "Hence, F_lower is updated as " << F_lower_updated << endl;
//    }
//    else if (break_method == 2){
//        cout << "Now, let us fix function limits to be divisible by beta. Originally F_lower: " << F_lower << " and F_upper: " << F_upper << endl;
//        current_length_upper = (F_upper - F_truth - 3*increments)/(0.5*(overline_k - 3)); //length of each piece in the upper half from the true value
//        desired_upper_frac = ceil(current_length_upper / increments); //desired length fraction of each piece
//        desired_upper = desired_upper_frac * increments; //desired length of each upper piece
//        cout << "The length of one upper interval is thus " << current_length_upper << " but desired length is: " << desired_upper << endl;
//        F_upper_updated = ((double) desired_upper*(0.5 * (overline_k - 3)) + (F_truth + 3*increments));
//        cout << "Hence, F_upper is updated as " << F_upper_updated << endl;
//        //below just to update lower neighboring limits
//        current_length_lower = (F_truth - increments - F_lower)/(0.5*(overline_k - 3)); //length of each piece in the upper half from the true value
//        desired_lower_frac = ceil(current_length_lower / increments); //desired length fraction of each piece
//        desired_lower = desired_lower_frac * increments; //desired length of each upper piece
//        cout << "The length of one lower interval is thus " << current_length_lower << " but desired length is: " << desired_lower << endl;
//        F_lower_updated = ((double) (F_truth - 3*increments) - desired_lower*(0.5 * (overline_k - 3)) );
//        cout << "Hence, F_lower is updated as " << F_lower_updated << endl;
//    }
//    else if (break_method == 3){
//        cout << "Now, let us fix function limits to be divisible by beta. Originally F_lower: " << F_lower << " and F_upper: " << F_upper << endl;
//        current_length_upper = (F_upper - F_truth - increments - f_overline)/(0.5*(overline_k - 5)); //length of each piece in the upper half from the true value
//        desired_upper_frac = ceil(current_length_upper / increments); //desired length fraction of each piece
//        desired_upper = desired_upper_frac * increments; //desired length of each upper piece
//        cout << "The length of one upper interval is thus " << current_length_upper << " but desired length is: " << desired_upper << endl;
//        F_upper_updated = ((double) desired_upper*(0.5 * (overline_k - 5)) + (F_truth + increments + f_overline));
//        cout << "Hence, F_upper is updated as " << F_upper_updated << endl;
//        //below just to update lower neighboring limits
//        current_length_lower = (F_truth - increments - f_overline - F_lower)/(0.5*(overline_k - 5)); //length of each piece in the upper half from the true value
//        desired_lower_frac = ceil(current_length_lower / increments); //desired length fraction of each piece
//        desired_lower = desired_lower_frac * increments; //desired length of each upper piece
//        cout << "The length of one lower interval is thus " << current_length_lower << " but desired length is: " << desired_lower << endl;
//        F_lower_updated = ((double) (F_truth - increments - f_overline) - desired_lower*(0.5 * (overline_k - 5)) );
//        cout << "Hence, F_lower is updated as " << F_lower_updated << endl;
//    }
//    else{
//        cout << "ISSUE with the data-dep lengths!" << endl;
//    }
    //increments = 0.02;
    double limits = 4.0;
    vector<double> X; arange(X, -limits, limits, increments); //input: construct the noise X -> here there is issue
    //
    int N = static_cast<int>(X.size()) - 1;
    double threshold = delta/1.20; //algorithmic hyperparameters (slack amount to delete)
    int option = 50, obj = 0, monoton = 1 ; //delete slack constraints every 500 iterations, objective is E[|X|]
    vector<vector<double>> p_decisions(overline_k, vector<double>(N)); //initialize the decision vectors (overline_k prob weight vectors)
    //below set the weights
    //uniform weiggt
    vector<double> obj_weight(overline_k);
    for (int k = 0; k < overline_k; k++) {
        obj_weight[k] = (double) 1.0/(1.0*(overline_k));
    }

    // middle four -- in case you would like to give more weights to the middle elements, etc
//    vector<double> obj_weight(overline_k);
//    for (int k = 0; k < overline_k; k++) {
//        if (k >= (overline_k/2 - 2) && k <= (overline_k/2 + 1)) {
//            obj_weight[k] = 0.25;
//        } else {
//            obj_weight[k] = 0;
//        }
//    }
//    print_vector(obj_weight);
    //triangular -- another weighting that gives more towards zero
//    vector<double> obj_weight(overline_k);
//    double middle = (overline_k - 1) / 2.0;
//    double weight_sum = 0.0;
//    for (int k = 0; k < overline_k; k++) {
//        double distance_from_middle = abs(k - middle);
//        obj_weight[k] = 1.0 - distance_from_middle / middle;
//        weight_sum += obj_weight[k];
//    }
//    for (int k = 0; k < overline_k; k++) {
//        obj_weight[k] /= weight_sum;
//    }
//
//    print_vector(obj_weight);
    //add more weight to the mid items
//    obj_weight[0] = 0.25;
//    obj_weight[obj_weight.size()-1] = 0.25;
//    obj_weight[obj_weight.size()-1] = 0.0;
//    obj_weight[(int) overline_k/2] +=1.0/(1.0*overline_k);
//    obj_weight[(int) overline_k/2 -1] +=1.0/(1.0*overline_k);
//    print_vector(obj_weight);
    //
    
    //    obj_weight[tr] = 1.0; // or in case only one is important
    //print benchmark results
    //double limit_of_trunc = (f_overline / epsilon) * log(1 + truncated_quantity);
    cout << "Truncated Laplace has support [-A, A] where A= " << (f_overline / epsilon) * log(1 + truncated_quantity) << endl;
    cout << "Truncated Laplace expected amplitude: " << (f_overline / epsilon) * (1 - (log(1 + truncated_quantity) / truncated_quantity)) << endl;
    cout << "Truncated Laplace expected power: " << 2 * pow((f_overline / epsilon), 2) * (1 - ((log(1 + truncated_quantity) + 0.5 * pow(log(1 + truncated_quantity), 2)) / truncated_quantity)) << endl;
    double duration = optimization_alg(X, F_lower_updated, F_upper_updated, f_overline, increments, overline_k, p_decisions, option, threshold, obj, break_method, obj_weight, monoton);
    cout << "Amplitude is " << noise_amplitude(X, obj_weight, p_decisions) << endl;
    cout << "Lower bound is " << noise_amplitude_lower(X, obj_weight, p_decisions) << endl;
    // push results
    vector<double> results;
    results.push_back(1.0*job_nr_to_read);
    results.push_back(f_overline);
    results.push_back(F_lower);
    results.push_back(F_upper);
    results.push_back(increments);
    results.push_back(1.0*overline_k);
    results.push_back(noise_amplitude(X, obj_weight, p_decisions));
    results.push_back(noise_amplitude_lower(X, obj_weight, p_decisions));
    results.push_back(duration/1000.0);
    //************************* Save Optimized Solution
    //two digit precision for the sensitivity of the function
//    std::stringstream stream; //start a stringstream
//    stream << std::fixed << std::setprecision(2) << f_overline; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
//    //two digit precision for the upper bound of the function
//    std::stringstream stream_upper; //start a stringstream
//    stream_upper << std::fixed << std::setprecision(2) << F_upper; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
//    //two digit precision for the lower bound of the function
//    std::stringstream stream_lower; //start a stringstream
//    stream_lower << std::fixed << std::setprecision(2) << F_lower; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
//
//    std::ofstream f; //start the file stream
//    f.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/new/X_sens" + stream.str() + "lower" + stream_lower.str() + "upper" + stream_upper.str()  + ".csv");
//    write_to_csv(X, f);
//    f.close();
//    //BELOW distirbution save
////    if (tr == 0){
////        f_dist.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/new/p_sens" + stream.str() + "lower" + stream_lower.str() + "upper" + stream_upper.str()  + ".csv");
////        write_to_csv(p_decisions[tr], f_dist); //note f_dist is constructed at the global scope
////    }
////    else{
////        write_to_csv(p_decisions[tr], f_dist);
////    }
////    if (tr == overline_k){
////    f_dist.close();
////    }
//    ABOVE distribution save
//    uncomment below
//    std::ofstream f_dist; //start the file stream
////    f_dist.open("/Users/.../Workspace/DP/Python/Visualise/multi/new_lhs.csv");
//    f_dist.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/SGD/dep/p_mult.csv");
//    for (int dtw = 0; dtw < overline_k; dtw++){
//        write_to_csv(p_decisions[dtw], f_dist); // made tr
//    }
////    write_to_csv(p_decisions[tr], f_dist); // made tr
//    f_dist.close();
//    f_dist.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/SGD/dep/X_mult.csv");
//    write_to_csv(X, f_dist); // made tr
//    f_dist.close();
//    //uncomment the above save
//    vector<double> underline_I(overline_k), overline_I(overline_k); //these will be filled soon]
//    if(break_method == 0){
//        domains(underline_I, overline_I, F_lower, F_upper_updated, overline_k, 0); // just to write the break-points to figure out which noise to take
//        f.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/new/ranges_sens" + stream.str() + "lower" + stream_lower.str() + "upper" + stream_upper.str()  + ".csv");
//        write_to_csv(overline_I, f);
//        f.close();
//    }
//    print_vector(p_decisions[2]);
    
    std::ofstream f;
    f.open("/Users/.../Workspace/DP/DPNBjulia/Results/Runtimes/mid_m_" + to_string(job_nr_to_read)+ ".csv");
    write_to_csv(results, f);
    f.close();
    
//    else{
//        domains_updated(underline_I, overline_I, F_lower_updated, F_upper_updated, F_truth, increments ,overline_k, break_method, f_overline); // just to write the break-points to figure out which noise to take
//        f.open("/Users/.../Workspace/DP/DPNBjulia/Distributions/new/ranges_sens" + stream.str() + "lower" + stream_lower.str() + "upper" + stream_upper.str()  + ".csv");
//        write_to_csv(overline_I, f);
//        f.close();
//    }
    return 0;
}

int main() {
    double f_overline = 2.0, F_lower = 0.0, F_upper = 4.0; //f_overline is the sensitivity Phi = [F_lower, F_upper]
//    int div_factor  = pow(2,0); // divide factor
    int job = 1; //1,2,3,4,5,6,7, ...
    overline_k = 4*pow(2,job-1); //int (F_upper - F_lower) - 2;//number of intervals (initially defined globally)
    double increments = 1/(1.0*pow(2,job-1)); //(F_upper-F_lower)/(1.0 *overline_k); //1.0/(1.0*div_factor); //increments = beta
    p.resize(overline_k); //just resize the decision variables
    main_multiple(f_overline, F_lower, F_upper, (int) (overline_k - 1)/2, increments, job); //tr points out which data-dep weight to take
    return 0;
}

// the following is a main function we used in HPC before. Keeping here for the sake of completeness to refer back later.
int main_hpc_old(int argc, const char * argv[], const char *param_env[]){
    //STEP 1: figure out the job number we are in
    int job_number = 0;
    for (int i = 0; param_env[i] != NULL; ++i){
        string str = param_env[i];
        if (str.find ("PBS_ARRAY_INDEX") != string::npos){
            job_number = atoi (str.substr (str.find_last_of ("=") + 1).c_str());
            break;
        }
    }
    double f_overline = 2.0, F_lower = 0.0, F_upper = 4.0;
    overline_k = 4*pow(2,job_number-1); //int (F_upper - F_lower) - 2;//number of intervals (initially defined globally)
    double increments = 1/pow(2,job_number-1);
    p.resize(overline_k); //just resize the decision variables
    main_multiple(f_overline, F_lower, F_upper, (int) (overline_k - 1)/2, increments, job_number); //tr points out which data-dep weight to take
    //STEP2: read the parameters array
//    std::string line;
//    std::vector<double> sensArray; //sensitivities
//    std::vector<double> minArray;  //F_lower's
//    std::vector<double> maxArray;  //F_upper's
//    // read senses
//    std::ifstream myfile("/rds/general/user/.../Datasets/sens_read_params.csv");
//    if(!myfile){
//        std::cout<<"Error opening output file"<< std::endl;
//        return -1;
//    }
//    while (std::getline(myfile, line, ',')){
//        sensArray.push_back(atof(line.c_str())); //push to the vector
//    }
//    // read mins
//    std::ifstream myfile_two("/rds/general/user/.../Datasets/mins_read_params.csv");
//    if(!myfile_two){
//        std::cout<<"Error opening output file"<< std::endl;
//        return -1;
//    }
//    while (std::getline(myfile_two, line, ',')){
//        minArray.push_back(atof(line.c_str())); //push to the vector
//    }
//    // read senses
//    std::ifstream myfile_three("/rds/general/user/.../Datasets/maxs_read_params.csv");
//    if(!myfile_three){
//        std::cout<<"Error opening output file"<< std::endl;
//        return -1;
//    }
//    while (std::getline(myfile_three, line, ',')){
//        maxArray.push_back(atof(line.c_str())); //push to the vector
//    }
    //STEP 3: now that we know the parameter, let us run that!
    //main_multiple(sensArray[job_number], minArray[job_number], maxArray[job_number]);
//    for (int tr = 1; tr < 2; tr++){
//        main_multiple(sensArray[job_number], minArray[job_number], maxArray[job_number], tr);//tr points out which data-dep weight to take
//        //now delete elemnts
//        for (int i = 0; i < p.size(); i++) {
//            p[i].clear();
//            //z[i].clear();
//        }
//    }
    //int div_factor = pow(2, job_number); //job nr: 0 -> 1, 1 -> 2, 2 -> 4, 3 -> 8, 4 -> 16, 5 -> 32, 6 -> 64, 7 -> 128.
//    double increments = 0.0; //dummy initialization //= 1.0/(1.0*div_factor);
//    //overline_k = 4*div_factor;
//
//    int job_counter = 0;
//    for (int i = 0; i < 8; i++){ //i is about the length of the 'F_k pieces'
//        if(job_counter == job_number){
//            int div_factor = pow(2,4);
//            increments = 1.0/(1.0*div_factor);
//            overline_k = 4*div_factor;
//            p.resize(overline_k);
//        }
//        job_counter++;
//    }
//    p.resize(overline_k);
//    main_multiple(2.0, 0.0, 4.0, (int) (overline_k - 1)/2, increments, job_number);
    
    return 0;
}
