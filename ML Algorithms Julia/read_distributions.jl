function read_once(epsilon, delta, f_overlines, F_lowers, F_uppers)
    Xs,ranges,ps,betas = Dict([]), Dict([]), Dict([]), Dict([]) #empty dicts
    for it in 1:length(f_overlines)
        lookup = sprintf1("%.2f", f_overlines[it])*"lower"*sprintf1("%.2f", F_lowers[it])*"upper"*sprintf1("%.2f", F_uppers[it])*".csv"
        X_temp = vec(readdlm("Distributions/new/X_sens"*lookup, ',', Float64))
        Xs[lookup] = X_temp
        betas[lookup] = round(X_temp[2] - X_temp[1],digits = 6) #take the beta
        ranges[lookup] = vec(readdlm("Distributions/new/ranges_sens"*lookup, ',', Float64))
        ps[lookup] = readdlm("Distributions/new/p_sens"*lookup, ',', Float64)
    end
    return Xs,ranges,ps,betas
end
