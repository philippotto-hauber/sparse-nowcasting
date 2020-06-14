#__________________________________________________________#
#_ Violin plots of absolute forecast error distribution 
#__________________________________________________________#


#____SET-UP________________________________________________#
Nh = 2

dir = "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/eval/sfe_distr/"

#______LOAD PACKAGES________________________________________#

using CSV
using Plots									# packages & plot backend		
using StatsPlots 
gr()

using Statistics							# some more packages
using StatsBase
using MLBase

#______PLOTS_______________________________________________#

function sort_df(df, N = 50_000)
    df_sort = df[1 : N, :] # this keeps the column names!
    for i in 1 : size(df, 2)
        tmp = sort(df[:, i])
        df_sort[:, i] = tmp[1 : N]
    end
    return df_sort
end 

function violin_plot(country, Nr, Nh)
    filename = "sfedistr_" * string(country) * "_Nr" * string(Nr) * "_Nh" * string(Nh)  * ".csv"

    df = CSV.read(dir * filename) # * concatenates strings

    # truncate draws
    #df = sort_df(df, 40_000)

    # square root of draws
    df = sqrt.(df)

    p = violin(["NIG"], df.BAR, alpha=0.8,side=:left,color=6,label="B-AR", title = "Nr = " * string(Nr), titlefontsize = 8)
    violin!(["PMNM"], df.BAR, alpha=0.8,side=:left,color=6,label="")
    violin!(["MG"], df.BAR, alpha=0.8,side=:left,color=6,label="")
    violin!(["HS"], df.BAR, alpha=0.8,side=:left,color=6,label="")
    
    violin!(["NIG"], df.NIG, alpha=0.8,side=:right,color=7,label="prior")
    violin!(["PMNM"], df.PMNM, alpha=0.8,side=:right,color=7,label="")
    violin!(["MG"], df.MG, alpha=0.8,side=:right,color=7,label="")
    violin!(["HS"], df.MG, alpha=0.8,side=:right,color=7,label="")

    return p
end

countries = ["GER", "US"]
for country in countries
    title = scatter(title = string(country) * ", recursive sample, surveys in levels", titlefontsize = 10, axis = nothing, border=:none, size = (25, 50))

    # plot 1
    Nr = 2
    p1 = violin_plot(country, Nr, Nh)
    # plot 2
    Nr = 5
    p2 = violin_plot(country, Nr, Nh)
    # plot 3
    Nr = 8
    p3 = violin_plot(country, Nr, Nh)

    plot(title, p1, p2, p3, layout = grid(4,1,heights=[0.001, 0.333, 0.333, 0.333]), size = (500, 1000), legend=:topleft, legendfontsize=4)

    # save figure
    savefig(dir * "fig_sfe_distr_" * country * ".pdf") 
end






