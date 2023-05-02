using Printf
function cpp(job)
    file_to_read = "./Results/Grid/figure1/figure_"*string(job)*".csv"
    return vec(readdlm(file_to_read, ','))
end

#check the grids first
#epsilon = 0
grids = zeros(100)
sub = zeros(100)
i = 1
#start with upper ones
for epsilon_row in 0:9
    job_nrs = (epsilon_row*10):((epsilon_row+1)*10 - 1)
    string_ub = ""
    string_lb = ""
    for job in job_nrs
        results = cpp(job)
        ub_achieved = results[4]    
        lb_achieved = results[6]
        o = (ub_achieved + lb_achieved)/2
        tlap_ub = results[2]
        tlap_lb = results[8]
        sub[i] = (ub_achieved - lb_achieved)/max(1,lb_achieved) #this is for our gaps
        print(tlap_ub)
        println(" "*string(o))
        val1 = (tlap_ub - o)/(max(o,1))
        # grids[i] = val1 #this is for our gaps
        val2 = (o-tlap_lb)/(max(o,1))
        grids[i] = val2
        i = i+ 1
    end
end
# max_grey = maximum(grids)
# min_grey = minimum(grids)

#epsilon = 0
epsilon_row = 9
job_nrs = (epsilon_row*10):((epsilon_row+1)*10 - 1)
string_ub = ""
string_lb = ""
for job in job_nrs
    results = cpp(job)
    ub_achieved = results[4]
    lb_achieved = results[6]
    o = (ub_achieved + lb_achieved)/2
    tlap_ub = results[2]
    println(tlap_ub)
    tlap_lb = results[8]
    val1 = (tlap_ub - o)/(max(o,1))
    val2 = ((o-tlap_lb)/max(o,1))
    grey_scale = mean(grids .<= val2)/(2.0)

    # (val1 + val2 - min_grey)/( (max_grey - min_grey) * 1.25) #normalize this guy to [0, 0.8] so that 1 - this one is in [0.2, 1]
    # string_ub = string_ub*"& "*"{\\cellcolor[gray]{"*string(round(1-grey_scale, digits = 3))*"} "*@sprintf("%.2f", (val1)*100)*"\\%}"
    string_lb = string_lb*"& "*"{\\cellcolor[gray]{"*string(round(1-grey_scale, digits = 3))*"} "*@sprintf("%.2f", val2*100)*"\\%}"
end
println(repeat("*",100))
println(string_ub)
println()
println(string_lb)



