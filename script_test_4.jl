
using CSV
using DataFrames
using JLD2

include("get_scenarios.jl")
include("optimization_one_comp.jl")
#utils 

MTBF = 1 # in years

nb_state = 4
P_evac = vcat([100.0 for i in 1:3], [0.0]) # vector of capacity levels

b = 3/(360*MTBF)

function build_P(theta)
    # theta = mean time in state 3

    P = zeros(4, 4, 2)

    for i in 1:4
        P[1,i,1] = 1.0 
    end

    mean_time_state_1 = (2/3)*360*MTBF - theta 
        
    P[1,1,2] = 1.0 - 1/mean_time_state_1
    P[2,1,2] = 1/mean_time_state_1
        
    P[3,3,2] = 1.0 - 1/theta
    P[4,3,2] = 1/theta

    P[2,2,2] = 1.0 - b
    P[3,2,2] = b

    P[4,4,2] = 1.0

    return P
end


# test

P_1 = build_P(365*MTBF/2)
P_2 = build_P(365*MTBF/100)


@assert all(x -> x <= 1, P_1)
@assert all(x -> 0 <= x, P_1)

tol = 1e-6  
@assert all(abs(sum(P_1[i, j, maint] for i in 1:nb_state) - 1) ≤ tol
            for j in 1:nb_state, maint in 1:2)

S = 5
years = 1
Q = 5
d = 5


Vals_all_1, Policies_m_1, Policies_k_1 = Bellman(years,Q,S,h,nb_state,P_1,d)
simulated_cost_11 = simulate(Policies_m_1,Policies_k_1,h,nb_state, years,Q,d,P_1)
simulated_cost_12 = simulate(Policies_m_1,Policies_k_1,h,nb_state, years,Q,d,P_2)

Vals_all_2, Policies_m_2, Policies_k_2 = Bellman(years,Q,S,h,nb_state,P_2,d)
simulated_cost_22 = simulate(Policies_m_2,Policies_k_2,h,nb_state, years,Q,d,P_1)
simulated_cost_21 = simulate(Policies_m_2,Policies_k_2,h,nb_state, years,Q,d,P_1)

filename_Vals_1 = "Policies/Vals_1.jld2"
@save filename_Vals_1 Vals_all_1
filename_Policies_m_1 = "Policies/Policies_m_1.jld2"
@save filename_Policies_m_1 Policies_m_1
filename_Policies_k_1 = "Policies/Policies_k_1.jld2"
@save filename_Policies_k_1 Policies_k_1

filename_Vals_2 = "Policies/Vals_2.jld2"
@save filename_Vals_2 Vals_all_2
filename_Policies_m_2 = "Policies/Policies_m_2.jld2"
@save filename_Policies_m_2 Policies_m_2
filename_Policies_k_2 = "Policies/Policies_k_2.jld2"
@save filename_Policies_k_2 Policies_k_2



# writing in CSV
header = DataFrame(simulated_cost_11=Float64[],simulated_cost_22=Float64[],simulated_cost_12=Float64[],simulated_cost_21=Float64[])
CSV.write("results_test_4.csv", header)
row = DataFrame(simulated_cost_11=simulated_cost_11,simulated_cost_22=simulated_cost_22,simulated_cost_12=simulated_cost_12,simulated_cost_21=simulated_cost_21)
CSV.write("results_test_4.csv", row; append=true, writeheader=false)
