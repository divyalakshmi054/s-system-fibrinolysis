include("Include.jl")
using DataFrames, Plots
using Plots.PlotMeasures
using StatsPlots
using CategoricalArrays

df = CSV.read(joinpath(pwd(),"sobol","Sensitivity-Sobol-1000-boot-100.csv"),DataFrame)
names = [ "α[1]"    ,  # 5
"α[2]"    ,  # 5
"α[3]"    ,  # 5
"α[4]"    ,  # 5
"α[5]"    ,  # 5
"G[FII, r1]"              ,  # 6
"G[FIIa, r1]"      ,  # 7
"G[FIIa, r2]"           ,  # 8
"G[AT, r2]"                 ,  # 9
"G[FIIa, 3]" ,
"G[FI, 3]"       ,  # 11
"G[tPA,4]"   , 
"G[Plgn,4]" ,
"G[PAI1,4]" ,
"G[Plasmin,5]"   ,
"G[TAFI,5]"   ,
"G[FIa,5]"];

total_order_i = df[:,:Total_order]
CI_total = df[:,:Total_order_CI]
first_order_i = df[:,:First_order]
CI_first = df[:,:First_order_CI]
comp_ind = vcat(total_order_i,first_order_i)
comp_CI = vcat(CI_total,CI_first)


_PATH_TO_FIGS = joinpath(pwd(),"figs")
_PATH_TO_SOBOL_FIGS = joinpath(_PATH_TO_FIGS,"sobol")
path_to_sobol_total_png = joinpath(_PATH_TO_SOBOL_FIGS, "total-sobol.png") 
path_to_sobol_total_pdf = joinpath(_PATH_TO_SOBOL_FIGS, "total-sobol.pdf") 
path_to_sobol_first_png = joinpath(_PATH_TO_SOBOL_FIGS, "first-sobol.png") 
path_to_sobol_first_pdf = joinpath(_PATH_TO_SOBOL_FIGS, "first-sobol.pdf") 
path_to_sobol_combo_png = joinpath(_PATH_TO_SOBOL_FIGS, "combo-sobol.png") 
path_to_sobol_combo_pdf = joinpath(_PATH_TO_SOBOL_FIGS, "combo-sobol.pdf") 

fig1 = bar(names, total_order_i, yerror=CI_total, yguidefont = 12, ms=8, msw = 1, c=:gray, bottom_margin = 20px, lw=1, xrotation=45, ylims = [0, 1.0], ylabel = "Total Order Sensitivity Index", legend = false, framestyle = :box)
savefig(fig1, path_to_sobol_total_png)
savefig(fig1, path_to_sobol_total_pdf)

fig2 = bar(names, first_order_i, yerror=CI_first, yguidefont = 12, ms=8, msw = 1, c=:gray, bottom_margin = 20px, lw=1, xrotation=45, ylims = [0, 1.0], ylabel = "First Order Sensitivity Index", legend = false, framestyle = :box)
savefig(fig2, path_to_sobol_first_png)
savefig(fig2, path_to_sobol_first_pdf)

function Base.unique(ctg::CategoricalArray)
    l = levels(ctg)
    newctg = CategoricalArray(l)
    levels!(newctg, l)
end

ctg = CategoricalArray(repeat(["Total Order", "First Order"], inner = 17))
levels!(ctg, ["Total Order", "First Order"])

nam = CategoricalArray(repeat(names, outer = 2))
levels!(nam, names)

fig3 = groupedbar(nam, comp_ind, yerror=comp_CI, group = ctg, yguidefont = 12, ms=5, msw = 1, c=[:gray :white], bottom_margin = 20px, lw=1, xrotation=45, ylims = [0, 1.0],ylabel = "Sensitivity Index", bar_width = 0.67, legend = true, framestyle = :box)
savefig(fig3, path_to_sobol_combo_png)
savefig(fig3, path_to_sobol_combo_pdf)


