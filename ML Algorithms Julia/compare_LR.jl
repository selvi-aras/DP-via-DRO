using JuMP #to call solvers and formulate optimization problems
using LinearAlgebra #to use functions like 'norm'
import Random #for random number generation
import MosekTools #MOSEK
import MathOptInterface #used by JuMP
const MOI = MathOptInterface #referring to MathOptInterface simply as 'MOI'
using MAT #if you want to read data from MATLABx
using Dualization
using DelimitedFiles #CSV, DataFrames,
using InvertedIndices
"Takes a model and adds a softplus constraint. See https://jump.dev/JuMP.jl/stable/tutorials/conic/logistic_regression"
function softplus(model, t, linear_transform) #exponential cone constraint
    z = @variable(model, [1:2], lower_bound = 0.0)
    #add the exp-cone constraints
    @constraint(model, sum(z) <= 1.0)
    @constraint(model, [linear_transform - t, 1, z[1]] in MOI.ExponentialCone())
    @constraint(model, [-t, 1, z[2]] in MOI.ExponentialCone())
end

"Returns a JuMP model to optimize for logistic regression."
function build_logit_model(X, y, regular, lambda)
    N, n = size(X) #N rows, n binary predictors
    model = Model(dual_optimizer(MosekTools.Optimizer)) #start the model via MOSEK
    set_optimizer_attributes(model, "MSK_IPAR_NUM_THREADS" => 1, "MSK_IPAR_INTPNT_MULTI_THREAD" => 0) #one thread
    @variable(model, beta[1:n]) #beta coefficients
    @variable(model, t[1:N]) #auxiliary variables
    for i in 1:N #add N softplus constraints, e.g., log-loss at i-th point <= t_i
        u = -(X[i, :]' * beta) * y[i]
        softplus(model, t[i], u)
    end
    # Define objective, which depends on whether we take regularization
    if regular == 1 #first order regularization
        @variable(model, 0.0 <= reg)
        @constraint(model, [reg; beta] in MOI.NormOneCone(n + 1))
        @objective(model, Min, sum(t)/N + (lambda * reg))
    else #no regularization
        @objective(model, Min, sum(t)/N)
    end
    return model
end
"Take data and return optimized logistic regression model."
function logistic_regression(X,y; regular = 0, lambda = 0)
    # Optimizes the logistic regression problem
    model = build_logit_model(X, y, regular, lambda)
    set_silent(model)
    JuMP.optimize!(model)
    if termination_status(model) != OPTIMAL #warn if not optimal
        error("Solution is not optimal.")
    end
    solver_time = solve_time(model)
    #see solution
    #beta_opt = JuMP.value.(model[:beta])
    #optimal_obj = JuMP.objective_value(model)
    return model, solver_time #question: does returning the model make things slow?
end


model, solver_time = logistic_regression(df_train, labels_train; regular = 1, lambda = 10^-8)
beta_lr = JuMP.value.(model[:beta])
in_lr, out_lr = mscs(beta_lr, df_train, df_test, labels_train, labels_test)
