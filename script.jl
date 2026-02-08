
using JLD2

include("get_scenarios.jl")
include("create_scenarios_prod.jl")

"""
---------------------------------------------------------------------------------------------
Parameters - can be modified 
---------------------------------------------------------------------------------------------
"""
nb_w = 500
nb_d = 1
years = 30
H = 6

"""
---------------------------------------------------------------------------------------------
Should not be modified
---------------------------------------------------------------------------------------------
"""

Q = 6

include("create_simple_pole.jl")
include("simulation multiple comps.jl")


#create_δ(C,h,pr,n,d,H,P_daily,"delta.jld2",nb_time,P_evac)
#create_δ(C,h_no_w,pr,n,d,H,P_daily,"delta_no_w.jld2",nb_time,P_evac)
#create_p(C,H,n,P_period,"p.jld2")

@load "delta99monopole.jld2"
#δ_w = δ 
#@load "delta_no_w.jld2"
#δ_no_w = δ 

@load "p99monopole.jld2"

policies_file = "policies_weakly2.csv"
open("policies_weakly2.csv", "w") do f
    # ouvrir en mode écriture tronque le fichier
end
policies_file_1 = "policies_benchmark_1.csv"
policies_file_2 = "policies_benchmark_2.csv"


simu_file = "simulation_weakly2.csv"
open("simulation_weakly2.csv", "w") do f
    # ouvrir en mode écriture tronque le fichier
end
simu_file_1 = "simulation_benchmark_1.csv"
simu_file_2 = "simulation_benchmark_2.csv"

cost_file = "cost_weakly2.csv"
open("cost_weakly2.csv", "w") do f
    # ouvrir en mode écriture tronque le fichier
end
cost_file_1 = "cost_benchmark_1.csv"
cost_file_2 = "cost_benchmark_2.csv"



# fluid approximation
cache = Dict{Tuple{Tuple{Vararg{Int64}}, Int64, Int64}, Tuple{Vector{Int64}, Int64}}()
collect_policies(policies_file,C)
simulate(nb_w,nb_d,H,h,pr,n,years,Q,d,P_daily,simu_file,policies_file,"fluid",δ,p,cost_file,false,C)

# do not delete otherwise keep policies in cache
#cache = Dict{Tuple{Tuple{Vararg{Int64}}, Int64, Int64}, Tuple{Vector{Int64}, Int64}}()

#heuristic (benchmark1)
#simulate(nb_w,nb_d,H,h,pr,n,years,Q,d,P_daily,simu_file_1,policies_file_1,"benchmark1",δ_w,p,cost_file_1,false,C)

# do not delete otherwise keep policies in cache
#cache = Dict{Tuple{Tuple{Vararg{Int64}}, Int64, Int64}, Tuple{Vector{Int64}, Int64}}()

# optimize without accessibility (benchmark2)
# collect_policies(policies_file_2, C)
# simulate(nb_w,nb_d,H,h,pr,n,years,Q,d,P_daily,simu_file_2,policies_file_2,"fluid",δ_no_w,p,cost_file_2,false,C)


#BuildPlot("costs.csv")
