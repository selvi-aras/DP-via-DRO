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

GRBEnv* env;        //GUROBI environment
GRBModel* model;    //model
vector<GRBVar> p;   //create a p vector, for probabilities, and the type of vector is gurobi var
GRBConstr* c;       //use later

double epsilon;     //will be fixed in the main fn
double delta;       //will be fixed in the main fn
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
            diff.push_back(v[i] - v[j]);
        }
    }
    return diff;
}

void arange(vector<double>& result, double start, double end, double increase) { //Numpy's arange copied, but limits included, so small modification
    result.clear();
    //result.reserve((end - start) / increase);
    for (double i = start; i <= end + pow(10,-7); i = i + increase) {
        result.push_back(round(i * pow(10,7)) / pow(10,7));
    }
}
//add a DP constraint that is given in the iterations of the algorithm
void add_dp_constr(GRBLinExpr& lhs, const vector<double>& X, int N, vector<double>& indexD, vector<double>& indexDD, vector<double>& intervals, double& index, double& indexd) {
    lhs = 0.0; //start building the lhs of constraint
    for (int i = 0; i < 2 * N + 1; i++) {
        index = indexD[i]; //index of i: [f(D) + x_i, f(D') + x_{i+1}] includes this interval
        if (index >= 0) {//else it means we cannot reach here, so don't add! //new addition: && index < p.size()-1 don't include to LHS if p(end) is in
            lhs += intervals[i] * (p[index] / (X[index + 1] - X[index]));//summation over i part of the constraint -> Len(A_ell)*p_i / Len(f(D) + x_i, f(D) + x_{i+1})
        }
        indexd = indexDD[i]; //index of j: [f(D') + x_j, f(D') + x_{j+1}] includes this interval
        if (indexd >= 0) {
            lhs -= exp(epsilon) * (intervals[i] * (p[indexd] / (X[indexd + 1] - X[indexd])));
        }
    }
    model->addConstr(lhs <= delta);
}
/*
**************************************
ESSENTIAL Functions Below (Algorithms)
**************************************
*/
//ALG3 below
void worst_event(const double epsilon, const double varphi, const vector<double>& X, const vector<GRBVar>& p, vector<double>& indexD, vector<double>& indexDD, vector<double>& intervals, double& violation) {
    //returns the worst-case violation for a fixed varphi (ALG3 of the paper!)
    //assume wlog that f(D) = 0 and f(D') = varphi
    double condition; //condition to check
    double length_ell; //length of each interval
    bool check_if_i; //checks if the interval of [Xj, X_j+1]  starts with an f(D)+i point or not
    violation = 0; //new addition (normally 0) //initially no violation
    int N = (int)X.size() - 1; // -1 because X has size N+1 and probability vector has size N (N intervals)
    int i = -1, j = -1; //position of the last noise index f(D) + x_i or f(D') + x_j. =-1 means no access or don't include this interval
    double f_under_ell; //this will keep the beginning of A_ell = [f_under_ell, f_over_ell], one interval from the common refinement
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
    double f_over_ell;
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
        //define length
        length_ell = f_over_ell - f_under_ell; //Len(A_ell)
        intervals[ell] = length_ell; //We will keep the interval lengths since these will be used in the privacy constraints (length * probability)
        //now check the three cases
            // case 1 (the first logical means that if the length of interval is zero don't add the probabilities since it will reduce sparsity of LP)
        if (abs(length_ell) <= pow(10, -7) || f_under_ell >= X[N] - pow(10, -7)  || f_over_ell <= X[0] + pow(10, -7)) { //new addition at the end
            // do nothing: don't add this interval to the worst-case A \subseteq Omega
            indexD[ell] = -1; //means no access
            indexDD[ell] = -1; //means no access
        }
        else if (f_under_ell >= X[N] + varphi - pow(10, -7) || f_over_ell <= X[0] + varphi + pow(10, -7)) { //means f(D') + X cannot access so always add!
            //add it!
            indexD[ell] = i;
            indexDD[ell] = -1;//-1 since it cannot access
            violation = violation + length_ell * (p[i].get(GRB_DoubleAttr_X) / (X[i + 1] - X[i])); //extra violation added by this A_\ell
        }
        else { //else check the condition
            condition = (p[i].get(GRB_DoubleAttr_X) / (X[i + 1] - X[i])) - exp(epsilon) * (p[j].get(GRB_DoubleAttr_X) / (X[j + 1] - X[j])); //condition to check
            if (condition > pow(10, -8)) { //then violation, so add this interval!
                indexD[ell] = i;
                indexDD[ell] = j;
                violation = violation + length_ell * condition;
            }
            else {
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
// meta algorithm below
double optimization_alg(const vector<double>& X, double f_overline, vector<double>& p_decisions, int option, double threshold, int most_viol, double amount, int obj, int monoton) {//decisions is a dummy to keep optimized p vector
    auto start = high_resolution_clock::now(); //start counting the time
    vector<double> Delta = differences(X); remove(Delta); upper_lower(Delta, f_overline); //set of potential worst-case [f(D') - f(D)] values, remove duplicates, remove the elements out of the range, add limits
    int N = (int)X.size() - 1;
    //dummies to be used later
    vector<double> indexD(2 * N + 1);    //unassigned -- all intervals
    vector<double> indexDD(2 * N + 1);   //unassigned
    vector<double> intervals(2 * N + 1); //unassigned
    double violation;
    double worst_violation;              //will keep worst violation out of all violations in Delta
    vector<double> worst_D(2 * N + 1);
    vector<double> worst_Dd(2 * N + 1);
    vector<double> worst_intervals(2 * N + 1);
    double index; double indexd;
    double slack;                       //will keep the slack for each constraints
    // ---- start GUROBI
    env = new GRBEnv(true); //start GUROBI environment
    //env->set("LogFile", "try_new_log.log"); //x-> operation is the same as *x.
    env->set(GRB_IntParam_OutputFlag, 0);   //mute GUROBI
    //env->set(GRB_IntParam_Method, 1);
    env->start();                           //I start here, before the while loop, which will make the model warm-start all the time!
    model = new GRBModel(*env);             // Start
    p.clear();                              //clear solution from the other iterations
    // Create variables
    for (int i = 0; i < N; i++) {//create p[0]... p[n-1]  (in total n)
        p.push_back(model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS)); //p is between 0-1
    }
    
    // Set objective: e.g., min amplitude (linear)!
    GRBLinExpr objchoice = 0;
    if (obj == 0) { //amplitude
        for (int i = 0; i < N; i++) {
            objchoice += p[i] * abs(X[i + 1] + X[i]) / 2;
        }
    }
    else if (obj == 1) { // power
        for (int i = 0; i < N; i++) {
            objchoice += p[i] * (pow(X[i + 1], 3) - pow(X[i], 3)) / (3 * (X[i + 1] - X[i]));
        }
    }
    else if (obj == 2) { // lower bound for amplitude
        for (int i = 0; i < N; i++) {
            objchoice += p[i] * min(abs(X[i]), abs(X[i + 1]));
        }
    }
    else if (obj == 3) { //an asymmetric distribution
        for (int i = 0; i < N/2; i++) {
            objchoice += p[i] * abs(X[i + 1] + X[i]) / 2;
        }
        for (int i = N/2; i <N ; i++) {
            objchoice += p[i] * abs(X[i + 1] + X[i]); //introduce an asymmetry
        }
    }
    
    if (obj == 4) { //new -- try to penalize the tails more and more
        for (int i = 0; i < N; i++) {
            objchoice += p[i] * pow(X[i+1], 10);
        }
    }
    model->setObjective(1 * objchoice, GRB_MINIMIZE);
    // Add constraint: sum(p) ==1 (most probably there is a shorter/faster implementation)
    GRBLinExpr sum_lhs = 0;
    for (int i = 0; i < N; i++) {
        sum_lhs += p[i];
    }
    model->addConstr(sum_lhs == 1); // probabilities sum up to 1
    if (monoton == 1) { //impose  monotonicity constraints
        for (int i = 0; i <= (N-2); i++) {
            if (X[i + 2] <= 0) {
                model->addConstr(p[i] / (X[i + 1] - X[i]) - p[i + 1] / (X[i + 2] - X[i + 1]) <= 0);
            }
            else if (X[i] >= 0) {
                model->addConstr(p[i] / (X[i + 1] - X[i]) - p[i + 1] / (X[i + 2] - X[i + 1]) >= 0);
            }
            else{
                model->addConstr(p[i] / (X[i + 1] - X[i]) - p[i + 1] / (X[i + 2] - X[i + 1]) == 0);
            }
        }
    }
//    else if (monoton == 0) { //add only symmetry-- do not comment out unless you want symmetry
//        for (int i = 0; i < floor(N/2); i++) {
//            model->addConstr(p[i] == p[N-i-1]);
//        }
//    }
    int iterations = 0; //each iteration of the while loop
    while (true) { //iterate algorithm until no more violation
        iterations++;
        worst_violation = 0; //start with zero (this variable will keep the worst violation out of all the worst-case event of each varphi \in Delta
        model->update();
        model->optimize();
        // slack deletion below |
        if (iterations % option == 0) {//the if condition ensures that slack deletion will be implemented every 'option'th steps
            c = model->getConstrs();
            for (int constiter = 0; constiter < model->get(GRB_IntAttr_NumConstrs); ++constiter) {
                slack = c[constiter].get(GRB_DoubleAttr_Slack);
                if (slack > threshold) {
                    model->remove(c[constiter]);
                }
            }
        }
        //check all \varphi \in \Delta below
        for (int varphi = 0; varphi < Delta.size(); varphi++) { //check worst case of all varphi \in \Delta
            worst_event(epsilon, Delta[varphi], X, p, indexD, indexDD, intervals, violation); //check the worst-case of the given varphi value (ALG3) -> update the violation parameter
            if (violation >= worst_violation + pow(10, -7)) { //if violation found in worst_event (ALG3) is largest so far, update!
                //NOTE: Here I am checking every solution and taking the WORST violated one.
                //save the worst-case scenario solution
                worst_D = indexD;
                worst_Dd = indexDD;
                worst_intervals = intervals;
                worst_violation = violation;
            }
            //below adds all violated constraints above a threshold if most_viol == 0 (activated)
            if (most_viol == 0 && violation > delta + amount - pow(10, -8)) { //then add all infaesible constraints
                add_dp_constr(sum_lhs, X, N, indexD, indexDD, intervals, index, indexd);
            }
        }
        //time to check whether the worst violation among all worst-case events is feasible or infeasible (stopping condition)
        if (worst_violation <= delta + pow(10,-8)) { //+0.0000001 is for solver precision
            break; //no more violations
        }
        else {//otherwise add the worst constraint
            add_dp_constr(sum_lhs, X, N, worst_D, worst_Dd, worst_intervals, index, indexd);
        }
    }
    auto stop = high_resolution_clock::now(); //end time
    auto duration = duration_cast<milliseconds>(stop - start); //time passed
    //****things to return
    for (int i = 0; i < N; i++) {//save all the GUROBI soltuions of p to the vector 'decisions'
        p_decisions[i] = p[i].get(GRB_DoubleAttr_X);
    }
    cout << "number of iterations: " << iterations << endl;
    cout << "time (ms): " << duration.count() << endl;
    return duration.count();
}

double noise_amplitude(vector<double> X, vector<double> p) { //in case you want to query
    double ampl = 0; //sum_i p_i \times (x_i + x_{i+1})/2
    for (int i = 0; i < p.size(); i++) {
        ampl += p[i] * abs(X[i] + X[i + 1]);
    }
    return (ampl/2);
}
double noise_amplitude_lower(vector<double> X, vector<double> p) {
    double ampl = 0; //sum_i p_i \times (x_i + x_{i+1})/2
    for (int i = 0; i < p.size(); i++) {
        ampl += p[i] * min(abs(X[i]) , abs(X[i + 1]));
    }
    return (ampl);
}
double noise_power(vector<double> X, vector<double> p) {
    double pp = 0;
    for (int i = 0; i < p.size(); i++) {
        pp += p[i] * (pow(X[i+1],3) - pow(X[i],3))/(X[i+1] - X[i]);
    }
    return (pp/3);
}
double noise_power_lower(vector<double> X, vector<double> p) {
    double pp = 0;
    for (int i = 0; i < p.size(); i++) {
        pp += p[i] * min(pow(X[i],2), pow(X[i+1],2));
    }
    return (pp);
}
/*
**************************************
Run and Test Functions
**************************************
*/
int main_single(double f_overline, double increments, int job_nr_to_read) { //tune and simulate a single experiment
    epsilon = 1.0; //desired epsilon
    delta =  0.1;   //desired delta
    //************************* Truncated Laplace
    double truncated_quantity = (exp(epsilon) - 1) / (2 * delta);

    double limits = 4.0; //our support L*beta is limits*beta
    //************************* Optimizing DP (comment below if you are giving increments manually/
    //Ignore the commented out codes -- something like this was used to
//    double limit_of_trunc = (f_overline / epsilon) *  log(1 + truncated_quantity);
//    double limits = ceil(limit_of_trunc); //noise varies from -limits to limits
//    cout << "limit taken: " << limits << endl;
//    double nr_vars = 1000.0;
//    double increments = 2*limits/nr_vars; //aka beta -- take it so that overall we have 1,000 decision variables
//    cout << "increments taken: " << increments << endl;
//    increments = f_overline/floor(f_overline/increments); //round the "increments" so that "f_overline" is divisible (we also state in the paper that Delta f / beta is divisible)
//    cout << "increments fixed for f_overline divisibility: " << increments << endl;
//    limits = nr_vars*increments / 2; //same rounding for [-Lim, +Lim] (in papers notation [-Lbeta, Lbeta]
//    cout << "Corrected limits: " << limits << endl;
    //
    vector<double> X; arange(X, -limits, limits, increments); //input: construct the noise X -> here there is issue
    //print_vector(X);
    //************************* What to Optimize
    int objchoice = 0; //Objective Function: 0: amplitude, 1: power, etc.
    int option = 100; //Constraint Deletion: in this many iterations we will delete most feasible constraints
    double threshold = delta/1.2; //Constraint Deletion: most feasible constraints will be deleted if they are more than this feasile
    int most_viol = 0; //Constrant Addition:  =0 means all of the violated constraints will be included, =1 means only the most violated one will be included
    double amount = delta/5.0; //Constrant Addition: this is the amount of violation where the corresponding constraint should be added if most_viol = 0
    int monoton = 0; //0: nothing, 1: impose monotonicity, etc.
    //************************* Truncated Stats
    //cout << "Gaussian Mechanism has amplitude = " << sqrt(2*log(1.25/delta)*pow(f_overline/epsilon, 2))*(sqrt(2)/sqrt(atan(1)*4)) << endl;
    vector<double> results;
    results.push_back((f_overline / epsilon) * log(1 + truncated_quantity));
    results.push_back((f_overline / epsilon) * (1 - (log(1 + truncated_quantity) / truncated_quantity)));
    results.push_back(2 * pow((f_overline / epsilon), 2) * (1 - ((log(1 + truncated_quantity) + 0.5 * pow(log(1 + truncated_quantity), 2)) / truncated_quantity)));
    cout << "Truncated Laplace has support [-A, A] where A= " << results[0] << endl;
    cout << "Truncated Laplace expected amplitude is: " << results[1] << endl;
    cout << "Truncated Laplace expected power is: " << results[2] << endl;
    double a = (delta + ((exp(epsilon) - 1)/2))/(exp(epsilon));
    double b = exp(-1.0*epsilon);
    //************************* Optimization Stats
    vector<double> p_decisions(X.size() - 1); //initializing the variable to be used in optimization_alg
    double duration = optimization_alg(X, f_overline, p_decisions, option, threshold, most_viol, amount, objchoice, monoton);
    results.push_back(noise_amplitude(X, p_decisions));
    results.push_back(noise_power(X, p_decisions));
    results.push_back(noise_amplitude_lower(X, p_decisions));
    cout << "Our amplitude is: " << results[3] << endl;
    cout << "Our power is: " << results[4] << endl;
    cout << "Amplitude lower is: " << results[5] << endl;
    double n = 0.0; //for next one
    //Find the TLAP LB
    double condition_sum = 0.0;
    while (condition_sum + a*pow(b,n) <= 0.50){//Theory in paper
        condition_sum += a*pow(b,(int)n);
        n++;
    }
    cout <<  "n is " << n << endl;
    results.push_back(2*a*f_overline*( (b - pow(b,n))/pow(1-b,2) - (n-1)*pow(b,n)/(1-b) )   );
    cout << "TLap Lower Bound: " << results[6] << endl;
    results.push_back(increments);
    cout << "All for increments: " << results[7] << endl;
    results.push_back(duration/1000.0);
    cout << "Duration (s): " << results[8] << endl;
    //cout << "Here is the optimal solution:" << endl;
    //print_vector(p_decisions);
    //************************* Save Optimized Solution
    std::stringstream stream; //start a stringstream
    stream << std::fixed << std::setprecision(2) << f_overline; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
    std::stringstream stream2; //start a stringstream
    stream2 << std::fixed << std::setprecision(1) << epsilon; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
    std::stringstream stream3; //start a stringstream
    stream3 << std::fixed << std::setprecision(1) << delta; // f_overline; // now take f_overline and convert it to 2 decimal string as we will make a file
    std::ofstream f;
//    f.open("/Users/.../Distributions/eps1.0delta0.1/X_use" + stream.str() + ".csv");
//    write_to_csv(X, f);
//    f.close();
//    f.open("/Users/.../Distributions/eps1.0delta0.1/p_use" + stream.str() + ".csv");
//    write_to_csv(p_decisions, f);
    //
//    std::ofstream f;
//    f.open("/Users/.../Results/Runtimes/Revised/right_" + to_string(job_nr_to_read)+ ".csv");
//    write_to_csv(results, f);
//    f.close();
    //

//    std::ofstream f;
//    f.open("/Users/.../DPNBjulia/Distributions/SGD/X_use_trying.csv"); // normally we dont have "trying" at the end, but 16 April 2023 I'm adding that to check better indep for SGD.
//    write_to_csv(X, f);
//    f.close();
//    f.open("/Users/.../DPNBjulia/Distributions/SGD/p_use_trying.csv");
//    write_to_csv(p_decisions, f);
//    f.close();
    return 0;
}

int main() {
    double f_overline = 1.0; //sensitivity of the query '\varphi' in paper
    int job_nr = 1; //job_nr -- just to save file
    double increments = 0.1; //'\beta' in paper
    main_single(f_overline, increments, job_nr); //call the data independent noise optimization function
    return 0;
}
// the following is the main function of HPC. Keeping here for the sake of completeness.
int main_hpc(int argc, const char * argv[], const char *param_env[]){
    //STEP 1: figure out the job number we are in
    double f_overline = 1.0;
    int job_number = 0;
    for (int i = 0; param_env[i] != NULL; ++i){
        string str = param_env[i];
        if (str.find ("PBS_ARRAY_INDEX") != string::npos){
            job_number = atoi (str.substr (str.find_last_of ("=") + 1).c_str());
            break;
        }
    }
    //STEP2: read the parameters array
    std::string line;
    std::vector<double> parArray;
    std::ifstream myfile("/rds/general/user/.../home/DPNBjulia/Datasets/betas_vector.csv"); //THIS IS TO READ betas and run the code for each
    if(!myfile){
        std::cout<<"Error opening output file"<< std::endl;
        return -1;
    }
    while (std::getline(myfile, line, ',')){
        parArray.push_back(atof(line.c_str())); //push to the vector
    }
    //STEP 3: now that we know the parameter, let us run that!
    main_single(f_overline,parArray[job_number],job_number);
    return 0;
}
