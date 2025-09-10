include("build_MC_converter.jl")
include("build_MC_cooling.jl")
include("utils_transition_probabilities.jl")
include("get_scenarios.jl")
include("create_scenarios_prod.jl")

using JLD2

#Parameters

Q = 6
T = 6
h_no_w =  ones(Int, 180, 1, 62)
h_optim = h # change for h_no_w to optimize without weather scenarios

t_start = 1
Tmax_computeδ = T+6

#components: converter, water circuit, fans, pumps, transformer, cooling, busbar
C = 7
n = [12, 3, 8, 3, 3, 3, 3]
d = [6, 1, 1, 2, 1, 1, 1] 

length_m = sort(unique(d))
insert!(length_m, 1, 0)
nb_time = length(length_m)

#initial state
x0 = [3 for c in 1:C]
k0 = Q

n_max = maximum(n) #we can create arrays with irregular size

P_evac = [100.0 for c in 1:C, x in 1:n_max]
for c in 1:C
    for x in n[c]:n_max
        P_evac[c,x] = 0.0
    end
end


function create_schedule_m(m,t_start,T_per)

    Policies_m_t = [0 for t in 1:T_per]
                    
    if m == 1
        Policies_m_t[t_start] = 1
    end

    return Policies_m_t

end

function create_schedule_κ(κ,t_start,T_per)

    Policies_κ_t = [0 for t in 1:T_per]

    if κ > 0
        for q in 1:κ
            Policies_κ_t[t_start + q - 1] = 1
        end
    end

    return Policies_κ_t

end

P_evac[3,1:8] = [100.0, 100.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0] #linear decrease for the fan

function Capacity(state_list,n,C,P_evac,q)
    if sum(q[c] for c in 1:C) >= 1 #at least one ongoing maintenance
        return 0
    else 
        return minimum([P_evac[c,state_list[c]] for c in 1:C])
    end
end

"""
---------------------------------------------------------------------------------------------
Transition probabilities
---------------------------------------------------------------------------------------------
"""

p = [0.0 for c in 1:C, x in 1:n_max, k in 0:Q, x′ in 1:n_max, k′ in 0:Q, m in 1:nb_time, κ in 0:Q, t in 1:Tmax_computeδ]

#transition probabilities for 2 months 


P_converter = build_P_converter(true_MTBF_converter,12)[:, :, 2]

MTBF_trans = 2222 
MTBF_cool = 122 
MTBF_bus = 4762 

P_trans = build_P_linear_daily(n[5], MTBF_trans)
P_cool = build_P_linear_daily(n[6], MTBF_cool)
P_bus = build_P_linear_daily(n[7], MTBF_bus)

P_daily = [P_converter,P_water,P_fan,P_pump,P_trans,P_cool,P_bus]
P_period = [P_converter^60,P_water^60,P_fan^60,P_pump^60,P_trans^60,P_cool^60,P_bus^60]

for c in 1:C
    for t in 1:Tmax_computeδ
        for k in 0:Q
            for x in 1:n[c]
                for κ in 0:k
                    for m in 1:nb_time
                        if length_m[m] >= d[c]
                            if (t % 6 == 0) #end of the year
                                p[c,x,k+1,1,Q + 1,m,κ + 1,t] = 1.0 #maintenance
                            else
                                p[c,x,k+1,1,k - κ + 1,m,κ + 1,t] = 1.0
                            end
                        else
                            for x′ in 1:n[c] 
                                if (t % 6 == 0) #end of the year
                                    p[c,x,k+1,x′,Q + 1,m,κ + 1,t] = P_period[c][x′,x] #no maintenance
                                else
                                    p[c,x,k+1,x′,k - κ + 1,m,κ + 1,t] = P_period[c][x′,x]
                                end
                            end
                        end
                    end

                end     
            end
        end   
    end
end

tol = 1e-10 #replace very small coefficients by 0 to help the solver

for i in eachindex(p)
    if abs(p[i]) <tol
        p[i] = 0.0
    end
end

@save "p.jld2" p

"""
---------------------------------------------------------------------------------------------
Costs
---------------------------------------------------------------------------------------------
"""

function average_period_cost(T,Policies_m_t,Policies_κ_t,state,h,pr,nb_state,d,P,P_evac)

    # simulate states only over a strategic period for one failure scenario, weather scenario given by h, starting from state

    u = [0 for i in 1:T]
    #ong_m = [0 for i in 1:T]
    q = [0 for i in 1:T]

    q[1] = d*round(Int,Policies_m_t[1])

    last_state = zeros(nb_state)
    last_state[state] = 1.0

    if (q[1]>0)&&(h[1]==1)
        u[1] = 1
    end

    costs = []

    if q[1] >= 1
        new_cost = pr[1]*(1-Policies_κ_t[1])
    else 
        new_cost = (1-Policies_κ_t[1])*max(pr[1] - sum(last_state[i]*P_evac[i] for i in 1:nb_state),0)
    end 
    push!(costs, new_cost)

    for t in 1:T-1

        if u[t] == 0
            new_state = P[:,:]*last_state # no maintenance
        else 
            new_state = [1.0; zeros(nb_state - 1)]  # maintenance
        end

        q[t+1] = q[t] - u[t] + d*round(Int,Policies_m_t[t+1])
        if (q[t+1]>0)&&(h[t+1]==1)
            u[t+1] = 1
        end

        if q[t+1] >= 1
            new_cost = pr[t]*(1-Policies_κ_t[t+1])
        else 
            new_cost =  (1-Policies_κ_t[t+1])*sum(max(pr[t] - P_evac[i], 0)*new_state[i] for i in 1:nb_state)
        end 

        last_state = new_state
        push!(costs, new_cost)
    end

    return sum(costs)

end 




δ = [0.0 for c in 1:C, x in 1:n_max, m in 1:nb_time, κ in 0:Q, t in 1:Tmax_computeδ]

nb_scen = size(h_optim, 2)

for c in 1:C
    for x in 1:n[c]
        for m_global in 1:nb_time
            for κ in 0:Q
                for t in 1:Tmax_computeδ

                    T_per = 62 - count(==( -1 ), h[t,1,:]) 

                    if length_m[m_global] >= d[c]
                        m = 1
                    else
                        m = 0
                    end

                    local Policies_m_t = create_schedule_m(m,t_start,T_per)
                    local Policies_κ_t = create_schedule_κ(κ,t_start,T_per)

                    local cost = 0.0
                    
                    for scen in 1:nb_scen #average on weather scenarios

                        local h_t = h_optim[t,scen,:]
                        local pr_t = pr[t,:]
                        cost += (1/nb_scen) * average_period_cost(T_per,Policies_m_t,Policies_κ_t,x,h_t,pr_t,n[c],d[c],P_daily[c],P_evac[c,:])

                    end

                    δ[c,x,m_global,κ+1,t] = cost

                end
            end
        end
    end
end

for i in eachindex(δ)
    if abs(δ[i]) <tol
        δ[i] = 0.0
    end
end


@save "delta.jld2" δ