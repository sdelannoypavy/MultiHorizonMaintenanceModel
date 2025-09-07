include("build_MC_converter.jl")
include("build_MC_cooling.jl")
include("utils_transition_probabilities.jl")
include("get_scenarios.jl")

#Parameters

Q = 20
T = 6

#components: converter, water circuit, fans, pumps, transformer, cooling, busbar
C = 7
n = [12, 3, 8, 3, 3, 3, 3]
d = [11, 1, 1, 2, 3, 1, 1] 

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

P_evac[3,1:8] = [100.0, 100.0, 100.0, 80.0, 60.0, 40.0, 20.0, 0.0] #linear decrease for the fan

"""
---------------------------------------------------------------------------------------------
Transition probabilities
---------------------------------------------------------------------------------------------
"""

p = [0.0 for c in 1:C, x in 1:n_max, k in 0:Q, x′ in 1:n_max, k′ in 0:Q, m in 0:1, κ in 0:Q, t in 1:T]

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
    for t in 1:T
        for k in 0:Q
            for x in 1:n[c]
                for κ in 0:k

                    if (t % 6 == 0) #end of the year
                        p[c,x,k+1,1,Q + 1,2,κ + 1,t] = 1.0 #maintenance
                    else
                        p[c,x,k+1,1,k - κ + 1,2,κ + 1,t] = 1.0
                    end

                    for x′ in 1:n[c] 
                        if (t % 6 == 0) #end of the year
                            p[c,x,k+1,x′,Q + 1,1,κ + 1,t] = P_period[c][x′,x] #no maintenance
                        else
                            p[c,x,k+1,x′,k - κ + 1,1,κ + 1,t] = P_period[c][x′,x]
                        end
                    end

                end     
            end
        end   
    end
end

"""
---------------------------------------------------------------------------------------------
Costs
---------------------------------------------------------------------------------------------
"""

function simulate_period(Policies_m_t,Policies_k_t,state,h,nb_state,d,P,P_evac)

    # simulate states only over a strategic period for one failure scenario, weather scenario given by h, starting from state

    T = 62 - count(==( -1 ), h)

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
        new_cost = 100*(1-Policies_k_t[1])
    else 
        new_cost = 100 - sum(last_state[i]*P_evac[i] for i in 1:nb_state)
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
            new_cost = 100*(1-Policies_k_t[t+1])
        else 
            new_cost = 100 - sum(new_state[i]*P_evac[i] for i in 1:nb_state)
        end 

        last_state = new_state
        push!(costs, new_cost)
    end

    if u[T] == 0
        last_state = P[:,:]*last_state # no maintenance
    else 
        last_state = [1.0; zeros(nb_state - 1)] # maintenance
    end

    return sum(costs)

end 




δ = [0.0 for c in 1:C, x in 1:n_max, m in 0:1, κ in 0:Q, t in 1:T]

nb_scen = size(h, 2)

for c in 1:C
    for x in 1:n[c]
        for m in 0:1
            for κ in 0:Q
                for t in 1:T

                    T_per = 62 - count(==( -1 ), h[t,1,:]) 

                    t_start = 30

                    Policies_m_t = [0 for t in 1:T_per]
                    Policies_m_t[t_start] = m

                    Policies_k_t = [0 for t in 1:T_per]

                    if κ > 0
                        for q in 1:κ
                            Policies_k_t[t_start + q - 1] = 1
                        end
                    end

                    local cost = 0.0
                    
                    for scen in 1:nb_scen #average on weather scenarios

                        h_t = h[t,scen,:]
                        cost += (1/nb_scen) * simulate_period(Policies_m_t,Policies_k_t,x,h_t,n[c],d[c],P_daily[c],P_evac[c,:])

                    end

                    δ[c,x,m+1,κ+1,t] = cost

                end
            end
        end
    end
end

