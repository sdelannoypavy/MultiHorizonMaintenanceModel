using Gurobi
using JuMP  
using Distributions
using Random
using StatsBase
#import MathOptInterface as MOI


# compute value function
function V(T,nb_maint,s,q_i,S,h,nb_state,P,Q,Vals,d,P_evac)

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    x_0 = zeros(nb_state)      
    x_0[s] = 1 

    @variable(model, 0 <= x[1:S,1:(T+1),1:nb_state] <= 1) # we need the value at T+1 to get next strategic period value function

    @variable(model, c[1:(T+1)] >= 0) # c[T+1] is the value function of the next strategic period

    @variable(model, m[1:T], Bin)
    @variable(model, m_type[1:nb_maint], Bin)
    @variable(model, u[1:S,1:T], Bin)
    @variable(model, ong_m[1:S,1:T], Bin) #ongoing maintenance
    @variable(model, q[1:S,1:T] >= 0)
    @variable(model, k[1:T], Bin)
    @variable(model, q_f, Int)
    @variable(model, δ[0:Q], Bin) # value of quota as the end of the period formulated using Bin. Usefull to write final cost.    
    @variable(model, P_m[1:nb_state, 1:nb_state] >=0)
    @variable(model, d_m >=0)

    @objective(model, Min, sum(c[t] for t in 1:(T+1)))

    M = 1e6  #Big M, could be changed to avoid numerical 
    alpha = 1 #multiplier for the code,could be increased to avoid numerical instability

    margin = 15 #margin to avoid maintenance that we can't finish
    @constraint(model, [t in 1:margin], m[T+1-t] == 0) 

    @constraint(model, sum(m_type[maint] for maint in 1:nb_maint) == 1)
    @constraint(model, sum(m[t] for t in 1:T) <= 1)
    @constraint(model, [i in 1:nb_state, j in 1:nb_state], P_m[i,j] == sum(m_type[maint]*P[i,j,maint] for maint in 1:nb_maint))

    # we consider that components are refirbushed since the first day of maintenance (no influence on cost), and no failure can happend during maintenance
    @constraint(model, [s in 1:S, t in 2:(T+1), i in 1:nb_state], x[s,t,i] == (1-ong_m[s,t-1])*sum(P[i,j,nb_maint]*x[s,t-1,j] for j in 1:nb_state) + ong_m[s,t-1]*sum(P_m[i,j]*x[s,t-1,j] for j in 1:nb_state))

    #@constraint(model, [s in 1:S, t in 1:(T)], sum(x[s,t,i] for i in 1:nb_state) == 1.0)
    @constraint(model, [s in 1:S, i in 1:nb_state], x[s,1,i] == x_0[i])

    @constraint(model, sum(m[t] for t in 1:T) <= 1) #maximum one maintenance for the strategic period 
    # PROBLEM: this contraint together with sum x = 1 leads to infeasibility (sometimes we want two maintenances a month...)

    # deterministic rule: "maintain as soon as possible"
    @constraint(model, d_m == sum(d[maint]*m_type[maint] for maint in 1:nb_maint))
    @constraint(model, [s in 1:S, i in 1:nb_state], q[s,1] == d_m*m[1])
    @constraint(model, [s in 1:S, t in 1:(T-1)], q[s,t+1] == q[s,t] - u[s,t] + d_m*m[t+1])

    @constraint(model, [s in 1:S, t in 1:T], u[s,t] + 1 - h[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], u[s,t] <= q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], u[s,t] <= h[s,t])

    @constraint(model, [s in 1:S, t in 1:T], ong_m[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], ong_m[s,t] <= q[s,t])

    @constraint(model, [t in 1:T], c[t] >= (alpha/S)*sum((1 - ong_m[s,t])*x[s,t,state]*(100 - P_evac[state]) + ong_m[s,t]*(1-k[t])*100*x[s,t,state] for s in 1:S, state in 1:nb_state)) # capacity equals 0 in state 12, maximum value everywhere else

    @constraint(model, q_f == q_i - sum(k[t] for t in 1:T)) 
    @constraint(model, q_f >= 0) 
    
    # ensure delta as the right value
    @constraint(model, [j in 0:Q], q_f - j <=  M * (1 - δ[j]))
    @constraint(model, [j in 0:Q], q_f - j >= -M * (1 - δ[j]))
    @constraint(model, sum(δ[j] for j in 0:Q) == 1)

    @constraint(model, [j in 0:Q], c[T+1] + M * (1 - δ[j]) >= (alpha/S)*sum(sum(x[s,T+1, i] * Vals[i, j+1] for i in 1:nb_state) for s in 1:S))

    optimize!(model)

    status = termination_status(model)
    # println("Statut de l'optimisation: $status")

    cost = objective_value(model)


    m = [value(m[t]) for t in 1:T]
    m_type = [value(m_type[maint]) for maint in 1:nb_maint]
    k = [value(k[t]) for t in 1:T]
    #u = [value(u[s,t]) for  s in 1:S, t in 1:T]
    q = [value(q[s,t]) for  s in 1:S, t in 1:T]
    ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    δ = [value(δ[j]) for j in 0:Q]

    #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    x = [value(x[1,t,state]) for  t in 1:(T+1), state in 1:nb_state]
    c = [value(c[t]) for t in 1:T]
    #sum_m = sum(value(m[t]) for t in 1:T)

    if status == MOI.OPTIMAL
        return(cost, m, m_type, k,q,δ,x,c,ong_m)
    else
        println("Aucune solution optimale trouvée.")
    end    

end


function Bellman(years,Q,S,h,nb_state,P)

    Tmax = 6*years

    Vals_old = zeros(nb_state, Q+1)

    # Remember all Bellman values
    Vals_all = Array{Float64, 3}(undef, nb_state, Q+1, Tmax)

    # Remember all optimal decisions
    Policies_m = Array{Float64}(undef, nb_state, Q+1, Tmax, 62)
    Policies_k = Array{Float64}(undef, nb_state, Q+1, Tmax, 62)

    for t in 1:Tmax
        Vals_new = zeros(nb_state, Q+1)
        for s in 1:nb_state
            for q in 0:Q
                h_t = h[Tmax-t+1, :, :]  # taille (S, T)
                T = 62 - count(==( -1 ), h_t[1,:]) # - 1 means ends of the month, so that we can represent months with variable lengths with vectors of the same dimensions
                val, m_opt, k_opt = V(T,nb_maint,s,q_i,S,h,nb_state,P,Q,Vals,d,P_evac)
                Vals_new[s, q+1] = val
                Policies_m[s,q+1,Tmax - t + 1,1:T] = m_opt
                Policies_k[s,q+1,Tmax - t + 1,1:T] = k_opt
            end
        end
        Vals_all[:, :, Tmax - t + 1] = Vals_new 
        Vals_old = Vals_new
    end

    return Vals_all, Policies_m, Policies_k

end

function simulate_period(Policies_m_t,state,h,nb_state)

    # simulate states only over a strategic period for one failure scenario, weather scenario given by h, starting from state

    states = [state]

    d = max(state - 1,1)

    T = 62 - count(==( -1 ), h)

    u = [0 for i in 1:T]
    q = [0 for i in 1:T]

    q[1] = d*Policies_m_t[1]

    last_state = state

    if (q[1]>0)&&(h[1]==1)
        u[1] = 1
    end

    for t in 1:T-1
        q[t+1] = q[t] - u[t] + d*Policies_m_t[t+1]
        if (q[t+1]>0)&&(h[t+1]==1)
            u[t+1] = 1
        end
        probas = [0.0 for i in 1:nb_state]

        if u[t] == 0
            probas = P_w[:, last_state]
        else 
            probas = P_m[:, last_state]
        end
        last_state = sample(1:nb_state, Weights(probas))

        push!(states, last_state)
    end

    if u[T] == 0
        probas = P_w[last_state, :]
    else 
        probas = P_m[last_state, :]
    end
    new_state = sample(1:nb_state, Weights(probas))

    return states,new_state

end 


function simulate(Policies_m,Policies_k,h,nb_state, years,Q)

    # simulate states over the entire concession period


    last_state = 12
    state_list = []
    last_q = Q

    for t in 1:(6*years)

        Policies_k_t = Policies_k[last_state,round(Int,last_q) + 1,t,:]
        sum_k = sum(Policies_k_t )

        Policies_m_t = Policies_m[last_state,round(Int,last_q) + 1,t,:]
        h_t = h[t, 1, :]

        new_states, last_state = simulate_period(Policies_m_t,last_state,h_t,nb_state)

        last_q -= sum_k
        
        append!(state_list, new_states)

    end

    return state_list

end


#script

#parameters

Q = 12 # quota of free maintenance days for each year 

S = 3

# random production scenario used for testing 
colonnes_neg1 = fill(-1, S, 2)  # 2 months of 30 days

#random wave height, accessible with proba 0.9
dist = Bernoulli(0.9)
h_month = [Int(rand(dist)) for s in 1:S, t in 1:60]
h_month = hcat(h_month, colonnes_neg1)

years = 1

#h = Array{Float64}(undef, 6*years, S, 62)
#for t in 1:(6*years)
#    h[t,:,:] = h_month
#end

#product = Array{Float64}(undef, 6*years, S, 62)
#for t in 1:(6*years)
#    product[t,:,:] = product_month
#end


#Vals_all, Policies_m, Policies_k = Bellman(years,Q,S,h,product,nb_state,P_w)
#state_list = simulate(Policies_m,Policies_k,h,nb_state, years,Q)

#@save "Vals_converter_120725.jld2" Vals_all
#@load "Vals_converter_120725.jld2" 
