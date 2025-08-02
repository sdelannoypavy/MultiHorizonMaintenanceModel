using Gurobi
using JuMP  
using Distributions
using Random
using StatsBase

#warning: this version of the code is specific to the converter

# compute value function
function V(T,s,q_i,S,h,product,nb_state,P_w,Q,Vals,is_end_year)

    x_0 = zeros(nb_state)      
    x_0[s] = 1 

    # maintenance duration depends on the state for the converter
    d = max(s - 1,1) 

    P_m = zeros(nb_state, nb_state)

    #P_m is defined here because it depends on d
    for i in 1:(d+1)
        P_m[1,i] = 1.0 
    end

    for i in (d+2):(nb_state) #we do not have enough maintenance days to maintain everything
        P_m[i-d,i] = 1.0 
    end

    model = Model(Gurobi.Optimizer)

    @variable(model, 0 <= x[1:S,1:(T+1),1:nb_state] <= 1) # we need the value at T+1 to get next strategic period value function

    @variable(model, c[1:(T+1)] >= 0) # c[T+1] is the value function of the next strategic period

    @variable(model, m[1:T], Bin)
    @variable(model, u[1:S,1:T], Bin)
    @variable(model, ong_m[1:S,1:T], Bin) #ongoing maintenance
    @variable(model, q[1:S,1:T] >= 0)
    @variable(model, k[1:T], Bin)
    @variable(model, q_f, Int)
    @variable(model, δ[0:Q], Bin) # value of quota as the end of the period formulated using Bin. Usefull to write final cost.    

    @objective(model, Min, sum(c[t] for t in 1:(T+1)))

    M = 1e6  #Big M, could be changed to avoid numerical 
    alpha = 1 #multiplier for the code,could be increased to avoid numerical instability

    margin = 5 #margin to avoid maintenance that we can't finish
    @constraint(model, [t in 1:(d + margin)], m[T+1-t] == 0) 

    # we consider that components are refirbushed since the first day of maintenance (no influence on cost), and no failure can happend during maintenance
    @constraint(model, [s in 1:S, t in 2:(T+1), i in 1:nb_state], x[s,t,i] == (1-ong_m[s,t-1])*sum(P_w[i,j]*x[s,t-1,j] for j in 1:nb_state) + ong_m[s,t-1]*sum(P_m[i,j]*x[s,t-1,j] for j in 1:nb_state))

    #@constraint(model, [s in 1:S, t in 1:(T)], sum(x[s,t,i] for i in 1:nb_state) == 1.0)
    @constraint(model, [s in 1:S, i in 1:nb_state], x[s,1,i] == x_0[i])

    @constraint(model, sum(m[t] for t in 1:T) <= 1) #maximum one maintenance for the strategic period 
    # PROBLEM: this contraint together with sum x = 1 leads to infeasibility (sometimes we want two maintenances a month...)

    # deterministic rule: "maintain as soon as possible"
    @constraint(model, [s in 1:S, i in 1:nb_state], q[s,1] == d*m[1])
    @constraint(model, [s in 1:S, t in 1:(T-1)], q[s,t+1] == q[s,t] - u[s,t] + d*m[t+1])

    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] + 1 - h[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] <= q[s,t])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] <= h[s,t])

    @constraint(model, [s in 1:S, t in 1:(T)], ong_m[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:(T)], ong_m[s,t] <= q[s,t])

    @constraint(model, [t in 1:T], c[t] >= (alpha/S)*sum((1 - ong_m[s,t])*x[s,t,nb_state]*product[s,t] + ong_m[s,t]*(1-k[t])*product[s,t] for s in 1:S)) # capacity equals 0 in state 12, maximum value everywhere else

    @constraint(model, q_f == q_i - sum(k[t] for t in 1:T)) 
    @constraint(model, q_f >= 0) 
    
    # ensure delta as the right value
    @constraint(model, [j in 0:Q], q_f - j <=  M * (1 - δ[j]))
    @constraint(model, [j in 0:Q], q_f - j >= -M * (1 - δ[j]))
    @constraint(model, sum(δ[j] for j in 0:Q) == 1)

    if is_end_year
        @constraint(model, [j in 0:Q], c[T+1] + M * (1 - δ[j]) >= (alpha/S)*sum(sum(x[s,T+1, i] * Vals[i, Q+1] for i in 1:nb_state) for s in 1:S)) #At year end, the quota resets
    else
        @constraint(model, [j in 0:Q], c[T+1] + M * (1 - δ[j]) >= (alpha/S)*sum(sum(x[s,T+1, i] * Vals[i, j+1] for i in 1:nb_state) for s in 1:S))
    end


    optimize!(model)

    status = termination_status(model)
    println("Statut de l'optimisation: $status")

    cost = objective_value(model)


    #numvars = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    #numcons = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))


    m = [value(m[t]) for t in 1:T]
    k = [value(k[t]) for t in 1:T]
    #u = [value(u[s,t]) for  s in 1:S, t in 1:T]
    q = [value(q[s,t]) for  s in 1:S, t in 1:T]
    δ = [value(δ[j]) for j in 0:Q]

    #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    #x = [value(x[s,t,nb_state]) for  s in 1:S, t in 1:(T+1)]
    #c = [value(c[t]) for t in 1:T]
    #sum_m = sum(value(m[t]) for t in 1:T)

    if status == MOI.OPTIMAL
        return(cost, m,k,q,δ)
    else
        println("Aucune solution optimale trouvée.")
    end    

end


function Bellman(years,Q,S,h,product,nb_state,P_w)

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
                product_t = product[Tmax-t+1, :, :]
                is_end_year = t % 6 == 1
                println(T,s,q)
                val, m_opt, k_opt = V(T,s,q,S,h_t,product_t,nb_state,P_w,Q,Vals_old,is_end_year)
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

p = 5.591277606933184e-5 # value for converter (see build_MC_converter.jl)
nb_state = 12

Q = 12 # quota of free maintenance days for each year 

S = 3

# random production scenario used for testing 
product_month = [100 for s in 1:S, t in 1:60]
colonnes_neg1 = fill(-1, S, 2)  # 2 months of 30 days
product_month = hcat(product_month, colonnes_neg1)

#random wave height, accessible with proba 0.9
d = Bernoulli(0.9)
h_month = [Int(rand(d)) for s in 1:S, t in 1:60]
h_month = hcat(h_month, colonnes_neg1)

P_w = zeros(nb_state, nb_state)

for i in 1:(nb_state - 1)
    P_w[i,i] = 1.0 - p 
    P_w[i+1,i] = p
end

P_w[nb_state,nb_state] = 1.0


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
