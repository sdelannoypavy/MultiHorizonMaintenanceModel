#to be used for tests with converter model

using CSV
using DataFrames
using JLD2
using Plots
using LaTeXStrings

include("get_scenarios.jl")
include("build_MC_converter.jl")
include("optimization_one_comp.jl")
include("get_params.jl")

header = DataFrame(MTBF_optim=Int[], cost=Float64[], simulated_cost=Float64[])
CSV.write("results_test_3.csv", header)

P_true = build_P_converter(true_MTBF)

for MTBF in MTBF_list 

    S = n_opt_list[1]
    years = years_list[1]
    Q = Q_list[1]
    d = d_list[1]

    local P = build_P_converter(MTBF)

    local Vals_all, Policies_m, Policies_k
    Vals_all, Policies_m, Policies_k = Bellman(years,Q,S,h,nb_state,P,d)
    local cost  = round(Int,Vals_all[1,Q+1,1])
    simulated_cost = round(Int,simulate(Policies_m,Policies_k,h,nb_state, years,Q,d,P_true))

    row = DataFrame(MTBF=MTBF, cost=cost, simulated_cost=simulated_cost)

    CSV.write("results_test_3.csv", row; append=true, writeheader=false)

    filename_Vals = "Policies/Vals_$(MTBF).jld2"
    @save filename_Vals Vals_all

    filename_Policies_m = "Policies/Policies_m_$(MTBF).jld2"
    @save filename_Policies_m Policies_m

    filename_Policies_k = "Policies/Policies_k_$(MTBF).jld2"
    @save filename_Policies_k Policies_k

end

println("✅ Results written in 'results_test_3.csv'")

# Chargement du CSV
df = CSV.read("results_test_3.csv", DataFrame)

plot(df.MTBF_optim, df.cost,
     xlabel="MTBF in optimization",
     ylabel="Simulated cost",
     legend=false,
     marker=:circle,
     label="Cost")

# Ajouter une barre verticale pointillée à true_MTBF
vline!([true_MTBF], linestyle=:dash, color=:red, label="true_MTBF")

annotate!(true_MTBF, -0.1, text(L"\theta^{\mathrm{MTBF}}", :red, 12, :center))


savefig("MTBF_sensitivity.pdf")
display(current())

