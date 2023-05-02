## Libs
#using PlotlyJS
using Measures
using Plots, LaTeXStrings
using Formatting, Printf
import GR
nr_to_sample = 5
## Privacy parameters
epsilon = 0.2
delta = 0.05
sensitivity = 2.0
var_Gauss = (2*log(1.25 / delta)*(sensitivity^2)) / (epsilon^2)
std_Gauss = sqrt(var_Gauss)
## values
ub_tlap = 4.72669 #from C++
lb_tlap = 3.91688 #from C++
ub_gauss = std_Gauss*sqrt(2)/sqrt(pi)
## try analytic gaussian
sigma = 6.653 #this is the value we need
quantity_1 = sensitivity/(2*sigma) - (epsilon*sigma)/(sensitivity)
quantity_2 = - sensitivity/(2*sigma) - (epsilon*sigma)/(sensitivity)
cdf(Normal(0.0,1.0), quantity_1) - exp(epsilon)*(cdf(Normal(0.0,1.0), quantity_2))
ub_analytic = sigma*sqrt(2)/sqrt(pi)
## reading our method
df = DataFrame(Support = Float64[], ubTLap = Float64[], powTLap = Float64[], ubOur = Float64[], powOur = Float64[], lbOur = Float64[], lbTlap = Float64[],increments = Float64[],runtime = Float64[])
for i in 1:nr_to_sample
    file_name = "./Results/Runtimes/Revised/right_"*string(i)*".csv"
    if isfile(file_name)
        push!(df, readdlm(file_name, ','))
    end
end
#store the "row number" in a column
df.row = [i for i in 1:size(df)[1]]
df.logruntime = log.(exp(1) .+ df.runtime)


df_mult = DataFrame(job = Float64[], sens = Float64[], Flower = Float64[], Fupper = Float64[], increments = Float64[], overlinek = Float64[], ubOur = Float64[],lbOur = Float64[], runtime = Float64[])
for i in 1:nr_to_sample
    file_name = "./Results/Runtimes/Revised/right_m_"*string(i)*".csv"
    if isfile(file_name)
        push!(df_mult, readdlm(file_name, ','))
    end
end
#store the "row number" in a column
df_mult.row = [i for i in 1:size(df_mult)[1]]
df_mult.logruntime = log.(exp(1) .+ df_mult.runtime)

#legend a little upper

plot_fig= Plots.plot((df.row), df.lbOur, fillrange = df.ubOur,  legend=:none,
        legend_font = 9, fmt = :pdf, xaxis=(L"1/\beta"*"; "*L"K/4", font(12)), xguidefontsize=50, yaxis=("expected "*L"\ell_1"*"-loss", font(12)),
        fillalpha = 0.35,xlims = (1,nr_to_sample), xticks=(1:nr_to_sample, [1,2,4,8,16]), ylims = (0.4,5.7), yticks = 0.5:1.0:5.5,
         c = 9, label = "Data Independent Noise", leftmargin = 4mm, rightmargin = 14mm, bottommargin = 5mm,topmargin = 5mm, dpi = 1000)
plot_fig= plot!((df.row),df.lbOur, msw = 0, ms = 1.5, c=9, label ="")
plot_fig= plot!((df.row),df.ubOur, msw = 0, ms = 1.5, c=9, label = "")

plot_fig = hline!([ub_analytic], alpha = 0.6, msw = 0, ms = 1.5,lw=2, c="black", label ="Analytic Gaussian")
plot_fig = hline!([ub_tlap], alpha = 0.6,msw = 0, ms = 1.5,lw=2, c="brown", label ="Truncated Laplace")
plot_fig = hline!([lb_tlap], alpha = 0.6,msw = 0, ms = 1.5,lw=2, c="magenta2", label = :"Near Optimal LB")
#now plot the runtime on the right-hand side
ax2 = twinx()
plot_fig = plot!(ax2, xguidefontsize=50, xlims = (1,nr_to_sample), xticks=-1:0, yaxis = ("runtime ("*L"\log(e + x)"*" seconds)", font(12)), linestyle = :dash, 
ylims = (0.4,7.2), yticks= 1:2:7, 
(df.row), df.logruntime, msw = 0, ms = 1.5, c = 9,  legend = :none, label = :none,leftmargin = 4mm, rightmargin = 14mm, bottommargin = 5mm,topmargin = 5mm, dpi = 1000)

# #multi
plot_fig= plot!((df_mult.row), df_mult.lbOur, fillrange = df_mult.ubOur,
        fillalpha = 0.35, c = 13, label = "Data Dependent Noise", dpi = 1000)
plot_fig= plot!((df_mult.row),df_mult.lbOur, msw = 0, ms = 1.5, c=13, label ="")
plot_fig= plot!((df_mult.row),df_mult.ubOur, msw = 0, ms = 1.5, c=13, label = "")
#now plot the runtime on the right-hand side
plot_fig = plot!(ax2, (df_mult.row),  df_mult.logruntime, linestyle = :dash, xlims = (1,nr_to_sample),xticks=-1:0,  msw = 0, ms = 1.5, c = 13,  legend = :none, label = :none)

# Plots.savefig(plot_fig, "Images/right_convergence_doublevertical.pdf")
df_mult.close = (df_mult.ubOur - df_mult.lbOur)./max.(1, (df_mult.ubOur - df_mult.lbOur))