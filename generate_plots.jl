include("Include.jl")
#pythonplot()
# which visit?
visit = 1
# tpa level?
tpa = 0

# loading up simulation data -
_PATH_TO_TMP_ = joinpath(pwd(),"tmp")
# load time vector from a sample simulation file (should all have same times)
time_df = CSV.read(joinpath(_PATH_TO_TMP_,"SIM-visit-$(visit)-Fib-$(tpa)-nM-tPA-run-1.csv"),DataFrame)
t = time_df[:,:Time]

sim_data = Array{Float64}(undef,(length(t),12))   # change number of columns to match number of sims/runs

for i ∈ 1:12
    _PATH_TO_FILE = joinpath(_PATH_TO_TMP_,"SIM-visit-$(visit)-Fib-$(tpa)-nM-tPA-run-$(i).csv")
    sim_df = CSV.read(_PATH_TO_FILE, DataFrame)
    temp = sim_df[:,:CF]
    sim_data[:,i] = temp
end


 _PATH_TO_REAL_DATA = joinpath(pwd(),"data")
real_df = CSV.read(joinpath(_PATH_TO_REAL_DATA,"Training-TEG-Transpose.csv"),DataFrame)
real_df = mapcols(col -> replace(col, missing => 0), real_df)
no_tpa_df = filter(:TPA => x->(x=="N"), real_df)
yes_tpa_df = filter(:TPA => x->(x=="Y"), real_df)

v1_no_df = filter(:Visit => x->(x=="V1"), no_tpa_df) 
v2_no_df = filter(:Visit => x->(x=="V2"), no_tpa_df) 
v4_no_df = filter(:Visit => x->(x=="V4"), no_tpa_df) 
v1_yes_df = filter(:Visit => x->(x=="V1"), yes_tpa_df) 
v2_yes_df = filter(:Visit => x->(x=="V2"), yes_tpa_df) 
v4_yes_df = filter(:Visit => x->(x=="V4"), yes_tpa_df) 

# dump sorted real data

real_t = collect(5:5:10800)/60
no_tpa_plot = transpose(Array{Float64}(no_tpa_df[:,4:end]))
yes_tpa_plot = transpose(Array{Float64}(yes_tpa_df[:,4:end]))
v1_no_tpa_plot = transpose(Array{Float64}(v1_no_df[:,4:end]))
v1_y_tpa_plot = transpose(Array{Float64}(v1_yes_df[:,4:end]))
v2_no_tpa_plot = transpose(Array{Float64}(v2_no_df[:,4:end]))
v2_y_tpa_plot = transpose(Array{Float64}(v2_yes_df[:,4:end]))
v4_no_tpa_plot = transpose(Array{Float64}(v4_no_df[:,4:end]))
v4_y_tpa_plot = transpose(Array{Float64}(v4_yes_df[:,4:end]))

mean_v1_no = mean(v1_no_tpa_plot,dims=2)
mean_v2_no = mean(v2_no_tpa_plot,dims=2)
mean_v4_no = mean(v4_no_tpa_plot,dims=2)

mean_v1_yes = mean(v1_y_tpa_plot,dims=2)
mean_v2_yes = mean(v2_y_tpa_plot,dims=2)
mean_v4_yes = mean(v4_y_tpa_plot,dims=2)

_PATH_TO_FIGS = joinpath(pwd(),"figs\\ens")
path_to_CF_ens_figs_png = joinpath(_PATH_TO_FIGS, "tPA_$(tpa)nM_visit_$(visit)_CF_ens.png") 
path_to_CF_ens_figs_pdf = joinpath(_PATH_TO_FIGS, "tPA_$(tpa)nM_visit_$(visit)_CF_ens.pdf") 

# plot set of simulations for a given visit number
fig1 = plot(t, sim_data, xticks=0.0:30:180, ylim = (0,80),label="",lw = 1,c=colorant"#89CCE2", bg="aliceblue", background_color_outside="white", framestyle = :box, xlabel="Time (min)", ylabel="CF (mm)", title="Clot firmness vs. time, visit $(visit), [tPA] = $(tpa)nM")
#    scatter!(real_t,v1_y_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")
    scatter!(real_t,v1_no_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")

#    scatter!(real_t,v2_y_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")
#    scatter!(real_t,v2_no_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")

#    scatter!(real_t,v4_y_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")
#    scatter!(real_t,v4_no_tpa_plot ,mc= :gray, msc= :gray, ma=0.2,ms=1,label="")

    savefig(fig1,path_to_CF_ens_figs_png)
    savefig(fig1,path_to_CF_ens_figs_pdf)







