
using JLD2

include("get_scenarios.jl")
include("parameters_weakly_coupled.jl")
years = 30
Q = 6

include("create_simple_pole.jl")
include("simulation multiple comps.jl")


function create_γ(C,h_optim,pr,n,d,T,P_daily,delta_file,nb_time,P_evac)
# γ remplaces δ

    n_max = maximum(n)

    γ = [0.0 for z in 1:2, m in 1:nb_time, κ in 0:Q, t in 1:T]

    nb_scen = size(h_optim, 2)

    for z in 1:2

        last_state = [1 for c in 1:C]

        if z == 2 #failure
            for c in 1:C
                if d[c] <= 2 #we suppose it is the maintenance with the lowest cost
                #if d[c] <= 2
                    last_state[c] = n[c]
                end
            end
        end

        for m_global in 1:nb_time
            for κ in 0:Q
                for t in 1:T

                    T_per = 62 - count(==( -1 ), h_optim[t,1,:]) 

                    local Policies_m_t_list = []

                    for c in 1:C

                        if length_m[m_global] >= d[c]
                            m = 1
                        else
                            m = 0
                        end

                        local Policies_m_t = create_schedule_m(m,t_start,T_per)
                        push!(Policies_m_t_list, Policies_m_t)
                    end

                    local Policies_κ_t = create_schedule_κ(κ,t_start,T_per)

                    local cost = 0.0
                    
                    for scen in 1:nb_scen #average on weather scenarios

                        local h_t = h_optim[t,scen,:]
                        local pr_t = pr[t,:]
                        cost_new, _ = simulate_period(Policies_m_t_list,Policies_κ_t,last_state,h_t,pr_t,n,d,P_daily)
                        cost += (1/nb_scen) * cost_new

                    end

                    γ[z,m_global,κ+1,t] = cost

                end
            end
        end

    end

    for i in eachindex(γ)
        if abs(γ[i]) <tol
            γ[i] = 0.0
        end
    end

    @save delta_file γ

end



function PLNE_borne_inf(H,C,x0,k0,t0,γ,p,n)

    Tmax = H

    α = create_α(x0,k0,C)

    model = Model(optimizer)
    #set_optimizer_attribute(model, "OutputFlag", 0)


    @variable(model, 0 <= q[c=1:C, x=1:n[c], k = 0:Q, m = 1:nb_time, κ = 0:Q, t = t0:Tmax] <= 1)     
    @variable(model, ρ[z=1:2, m = 1:nb_time, κ = 0:Q, t = t0:Tmax] >= 0)    
    @variable(model, A[κ = 0:Q, m = 1:nb_time, t = t0:Tmax] >= 0)   
    @variable(model, μ[k = 0:Q, t = t0:Tmax] >= 0)  

    M = 1.0

    @objective(model, Min, M*sum(γ[z,m,κ+1,t]*ρ[z,m,κ,t] for z in 1:2, m in 1:nb_time, κ in 0:Q, t in t0:Tmax))

    # do we really want one vector of probability transitions for each t?
    @constraint(model, [c = 1:C, t = (t0+1):Tmax, x = 1:n[c], k = 0:Q], sum(q[c,x,k,m,κ,t] for m in 1:nb_time, κ in 0:Q) 
        >= sum(p[c,x′,k′+1,x,k+1,m,κ+1,t]*q[c,x′,k′,m,κ,t-1] for x′ in 1:n[c], k′ in 0:Q, m in 1:nb_time, κ in 0:Q))

    @constraint(model, [c = 1:C, m = 1:nb_time, κ = 0:Q, t = t0:Tmax], sum(q[c,x,k,m,κ,t] for x in 1:n[c], k in 0:Q) == A[κ,m,t])

    @constraint(model, [c = 1:C, x = 1:n[c], k = 0:Q], sum(q[c,x,k,m,κ,t0] for m in 1:nb_time, κ in 0:Q) == α[c,x,k+1])

    @constraint(model, [c = 1:C, x = 1:n[c], k = 0:(Q-1), m = 1:nb_time, κ = (k+1):Q, t = t0:Tmax], q[c,x,k,m,κ,t] == 0)

    @constraint(model, [c = 1:C, k = 0:Q, t = t0:Tmax], sum(q[c,x,k,m,κ,t] for x in 1:n[c], m in 1:nb_time, κ in 0:Q) == μ[k,t])

    @constraint(model, [c = 1:C, m = 1:nb_time, κ = 0:Q, t = t0:Tmax], sum(q[c,n[c],k,m,κ,t] for k in 0:Q) <= ρ[2,m,κ,t])
    @constraint(model, [m = 1:nb_time, κ = 0:Q, t = t0:Tmax],  ρ[1,m,κ,t] == A[κ,m,t] - ρ[2,m,κ,t])


    """
    ---------------------------------------------------------------------------------------------
    Optimize
    ---------------------------------------------------------------------------------------------
    """

    optimize!(model)

    status = termination_status(model)
    cost = objective_value(model)

    q_val = [0.0 for c in 1:C, x in 1:maximum(n), k in 0:Q, m in 1:nb_time, κ in 0:Q, t in 1:(Tmax - t0 + 1)]

    for c in 1:C
        for x in 1:n[c], k in 0:Q, m in 1:nb_time, κ in 0:Q, t in t0:Tmax
            q_val[c, x, k+1, m, κ+1, t - t0 + 1] = value(q[c, x, k, m, κ, t])
        end
    end

    return cost, q_val

end

k0 = Q
x0 = [1 for c in 1:C]

T = years*6

#create_p(C,T,n,P_period,"p_30.jld2") #years-6 because of preprocessing in the function 
#create_γ(C,h_optim,pr,n,d,T,P_daily,"delta_30_inf.jld2",nb_time,P_evac)

@load "p_30.jld2"
@load "delta_30_inf.jld2"

cost,q = PLNE_borne_inf(T,C,x0,k0,1,γ,p,n)