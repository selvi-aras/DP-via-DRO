using SpecialFunctions

function Phi(t)
    return 0.5*(1.0 + erf(float(t)/sqrt(2.0)))
end

function caseA(epsilon,s)
    return Phi(sqrt(epsilon*s)) - exp(epsilon)*Phi(-sqrt(epsilon*(s+2.0)))
end

function caseB(epsilon,s)
    return Phi(-sqrt(epsilon*s)) - exp(epsilon)*Phi(-sqrt(epsilon*(s+2.0)))
end

function doubling_trick(predicate_stop, s_inf, s_sup)
    while !(predicate_stop(s_sup))
        s_inf = s_sup
        s_sup = 2.0*s_inf
    end
    return s_inf, s_sup
end

function binary_search(predicate_stop, predicate_left, s_inf, s_sup)
    s_mid = s_inf + (s_sup-s_inf)/2.0
    while !(predicate_stop(s_mid))
        if (predicate_left(s_mid))
            s_sup = s_mid
        else
            s_inf = s_mid
        end
        s_mid = s_inf + (s_sup-s_inf)/2.0
    end
    return s_mid
end

"Julia version of https://github.com/BorjaBalle/analytic-gaussian-mechanism/blob/master/agm-example.py"
function calibrateAnalyticGaussianMechanism(epsilon, delta, GS; tol = 10^(-12))
    delta_thr = caseA(epsilon, 0.0)

    if (delta == delta_thr)
        alpha = 1.0
    else
        if (delta > delta_thr)
            predicate_stop_DT = s -> caseA(epsilon, s) >= delta
            function_s_to_delta = s-> caseA(epsilon, s)
            predicate_left_BS = s -> function_s_to_delta(s) > delta
            function_s_to_alpha = s -> sqrt(1.0 + s/2.0) - sqrt(s/2.0)
        else
            predicate_stop_DT = s ->  caseB(epsilon, s) <= delta
            function_s_to_delta = s ->  caseB(epsilon, s)
            predicate_left_BS = s ->  function_s_to_delta(s) < delta
            function_s_to_alpha = s ->  sqrt(1.0 + s/2.0) + sqrt(s/2.0)
        end

        predicate_stop_BS = s -> abs(function_s_to_delta(s) - delta) <= tol

        s_inf, s_sup = doubling_trick(predicate_stop_DT, 0.0, 1.0)
        s_final = binary_search(predicate_stop_BS, predicate_left_BS, s_inf, s_sup)
        alpha = function_s_to_alpha(s_final)
    end
    sigma = alpha*GS/sqrt(2.0*epsilon)

    return sigma
end

# "Sample noise from the analytic Gaussian mechanism."
# function analytic_noise(sigma_Gauss::Float64)
#     if sigma_Gauss < 0
#         error("Revise the parameters.")
#     end
#     return rand(Normal(0.0, sigma_Gauss), 1)[1] #sample from normal dist
# end
