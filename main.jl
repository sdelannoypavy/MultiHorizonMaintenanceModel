#to be used for tests with converter model

using CSV
using DataFrames
using JLD2

include("get_scenarios.jl")
include("build_MC_converter.jl")
include("optimization_one_comp.jl")
include("get_params.jl")

header = DataFrame(MTBF=Int[], S=Int[], years=Int[], Q=Int[], d=Int[], cost=Float64[], simulated_cost=Float64[])
CSV.write("results.csv", header)

for MTBF in MTBF_list 

    if !(MTBF === nothing) 
        global P = build_P_converter(MTBF)
    else
        MTBF = true_MTBF #to get the right value in the CSV file
    end

    for d in d_list
        for S in n_opt_list
            for years in years_list
                for Q in Q_list

                    for n_test in 1:nb_test_per_data

                        local Vals_all, Policies_m, Policies_k
                        Vals_all, Policies_m, Policies_k = Bellman(years,Q,S,h,nb_state,P,d)
                        local cost  = Vals_all[1,Q+1,1]
                        simulated_cost = simulate(Policies_m,Policies_k,h,nb_state, years,Q,d,P)

                        row = DataFrame(MTBF=MTBF, S=S, years=years, Q=Q, d=d, cost=cost, simulated_cost=simulated_cost)

                        CSV.write("results.csv", row; append=true, writeheader=false)

                        filename_Vals = "Policies/Vals_$(d)_$(S)_$(years)_$(Q)_$(n_test)_$(MTBF).jld2"
                        @save filename_Vals Vals_all

                        filename_Policies_m = "Policies/Policies_m_$(d)_$(S)_$(years)_$(Q)_$(n_test)_$(MTBF).jld2"
                        @save filename_Policies_m Policies_m

                        filename_Policies_k = "Policies/Policies_k_$(d)_$(S)_$(years)_$(Q)_$(n_test)_$(MTBF).jld2"
                        @save filename_Policies_k Policies_k

                    end

                end
            end
        end
    end
end

println("✅ Results written in 'results.csv'")

