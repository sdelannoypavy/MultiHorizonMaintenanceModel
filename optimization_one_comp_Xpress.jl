using Xpress
using JuMP
using Distributions
using Random
using StatsBase

Xpress.Lib.XPRSinit("")


# compute value function
function V(T,s_init,q_i,S,h,nb_state,P,Q,Vals,d,P_evac)

    println(s_init,q_i)
    model = Model(Xpress.Optimizer)

    x_0 = zeros(nb_state)      
    x_0[s_init] = 1 

    tol = 1e-3

    @variable(model, 0 <= x[1:S,1:(T+1),1:nb_state] <= 1 + tol) # we need the value at T+1 to get next strategic period value function

    @variable(model, c[1:(T+1),1:nb_state] >= 0) # c[T+1] is the value function of the next strategic period

    @variable(model, m[1:T], Bin)
    @variable(model, u[1:S,1:T], Bin)
    @variable(model, ong_m[1:S,1:T], Bin) #ongoing maintenance
    @variable(model, q[1:S,1:T] >= 0)
    @variable(model, k[1:T], Bin)
    @variable(model, q_f, Int)
    @variable(model, δ[0:Q], Bin) # value of quota as the end of the period formulated using Bin. Usefull to write final cost. 

    M = 1e6  #Big M, could be changed to avoid numerical 
    alpha = 1.0

    @objective(model, Min, (alpha/S)*sum(c[t,s] for t in 1:(T+1), s in 1:nb_state))

    margin = 30 #margin to avoid maintenance that we can't finish
    @constraint(model, [t in 1:margin], m[T+1-t] == 0) 
 
    # we consider that components are refirbushed since the first day of maintenance (no influence on cost), and no failure can happend during maintenance
    @constraint(model, [s in 1:S, t in 2:(T+1), i in 1:nb_state], !ong_m[s,t-1] --> {x[s,t,i] == sum(P[i,j,2]*x[s,t-1,j] for j in 1:nb_state)})
    @constraint(model, [s in 1:S, t in 2:(T+1), i in 1:nb_state], ong_m[s,t-1] --> {x[s,t,i] == sum(P[i,j,1]*x[s,t-1,j] for j in 1:nb_state)})

    #@constraint(model, [s in 1:S, t in 1:(T)], sum(x[s,t,i] for i in 1:nb_state) == 1.0)
    @constraint(model, [s in 1:S, i in 1:nb_state], x[s,1,i] == x_0[i])
 

    #@constraint(model, [s in 1:S, t in 1:(T+1)], sum(x[s,t,state] for state in 1:nb_state) == 1)


    # deterministic rule: "maintain as soon as possible"

    @constraint(model, [s in 1:S, i in 1:nb_state], q[s,1] == d*m[1])
    @constraint(model, [s in 1:S, t in 1:(T-1)], q[s,t+1] == q[s,t] - u[s,t] + d*m[t+1])

    @constraint(model, [s in 1:S, t in 1:T], u[s,t] + 1 - h[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], u[s,t] <= q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], u[s,t] <= h[s,t])

    @constraint(model, [s in 1:S, t in 1:T], ong_m[s,t] >= (1/T)*q[s,t])
    @constraint(model, [s in 1:S, t in 1:T], ong_m[s,t] <= q[s,t])

    @constraint(model, [s in 1:S, t in 1:T], !ong_m[s,t] --> {c[t,s] == sum(x[s,t,state]*(100 - P_evac[state]) for state in 1:nb_state)})
    @constraint(model, [s in 1:S, t in 1:T], ong_m[s,t] --> {c[t,s] == (1-k[t])*100})


    @constraint(model, q_f == q_i - sum(k[t] for t in 1:T)) 
    @constraint(model, q_f >= 0)

    # ensure delta as the right value
    @constraint(model, [j in 0:Q], q_f - j <=  M * (1 - δ[j]))
    @constraint(model, [j in 0:Q], q_f - j >= -M * (1 - δ[j]))
    @constraint(model, sum(δ[j] for j in 0:Q) == 1)
 
    
    @constraint(model, [j in 0:Q, s in 1:S], δ[j] --> {c[T+1,s] == sum(x[s,T+1, i] * Vals[i, j+1] for i in 1:nb_state)})


    optimize!(model)

    status = termination_status(model)
    # println("Statut de l'optimisation: $status")

    cost = objective_value(model)


    m = [value(m[t]) for t in 1:T]
    k = [value(k[t]) for t in 1:T]
    #u = [value(u[s,t]) for  s in 1:S, t in 1:T]
    #q = [value(q[s,t]) for  s in 1:S, t in 1:T]
    #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    #δ = [value(δ[j]) for j in 0:Q]

    #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    #x = [value(x[1,t,state]) for  t in 1:(T+1), state in 1:nb_state]
    #c = [value(c[t]) for t in 1:T]
    #sum_m = sum(value(m[t]) for t in 1:T)

    if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
        return(cost, m, k)
    else
        println("Aucune solution optimale trouvée.")
        println(status)
    end    

end


function Bellman(years,Q,S,h,nb_state,P,d)

    Tmax = 6*years

    h_reduced = reduce_array(h,S)#randomly select S scenarios

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
                h_t = h_reduced[Tmax-t+1, :, :]  # taille (S, T)
                T = 62 - count(==( -1 ), h_t[1,:]) # - 1 means ends of the month, so that we can represent months with variable lengths with vectors of the same dimensions
                println(t)
                val, m_opt, k_opt = V(T,s,q,S,h_t,nb_state,P,Q,Vals_old,d,P_evac)
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

function simulate_period(Policies_m_t,Policies_k_t,state,h,nb_state,d)

    # simulate states only over a strategic period for one failure scenario, weather scenario given by h, starting from state

    states = [state]

    T = 62 - count(==( -1 ), h)

    u = [0 for i in 1:T]
    q = [0 for i in 1:T]

    q[1] = d*Policies_m_t[1]

    new_state = state

    if (q[1]>0)&&(h[1]==1)
        u[1] = 1
    end

    costs = []
    new_cost = 0.0

    for t in 1:T-1
        q[t+1] = q[t] - u[t] + d*Policies_m_t[t+1]
        if (q[t+1]>0)&&(h[t+1]==1)
            u[t+1] = 1
        end
        probas = [0.0 for i in 1:nb_state]

        if u[t] == 0
            probas = P[:, new_state,1]
        else 
            probas = P[:, new_state, 2]
        end
        new_state = sample(1:nb_state, Weights(probas))

        if q[t] >= 1
            k[t] = Policies_k_t[t]
            new_cost = 100*(1-k[t])
        else 
            new_cost = 100 - P_evac[new_state]
        end 

        push!(states, new_state)
        push!(costs, new_cost)
    end

    if u[T] == 0
        probas = P[new_state, :, 1]
    else 
        probas = P[new_state, :, 2]
    end
    new_state = sample(1:nb_state, Weights(probas))

    return sum(costs), states,new_state

end 


function simulate_concession_period(Policies_m, Policies_k, h, nb_state, years, Q, d)

    # simulate states over the entire concession period

    total_cost = 0

    last_state = 12
    last_q = Q

    for t in 1:(6*years)

        Policies_k_t = Policies_k[last_state,round(Int,last_q) + 1,t,:]
        sum_k = sum(Policies_k_t )

        Policies_m_t = Policies_m[last_state,round(Int,last_q) + 1,t,:]
        h_t = h[t, :]

        cost, new_states, last_state = simulate_period(Policies_m_t, Policies_k_t, last_state, h_t, nb_state, d)

        last_q -= sum_k
        
        total_cost += cost

    end

    return total_cost

end



function simulate(Policies_m,Policies_k,h,nb_state, years,Q, d)

    list_cost = []

    nb_scen = size(h, 2)

    for scen in 1:nb_scen 
        new_cost = simulate_concession_period(Policies_m, Policies_k, h[:,scen,:], nb_state, years, Q, d)
        push!(list_cost, new_cost)
    end

    return mean(list_cost)

end


function reduce_array(A::Array{T,3}, l::Int) where T
    n, m, k = size(A)
    @assert l ≤ m "l doit être inférieur ou égal à m"

    # Tirer au hasard l indices uniques parmi 1:m
    selected_cols = sort(sample(1:m, l; replace=false))

    # Extraire les colonnes sélectionnées
    A_reduced = A[:, selected_cols, :]

    return A_reduced
end



#Vals_all, Policies_m, Policies_k = Bellman(years,Q,S,h,nb_state,P)
#state_list = simulate(Policies_m,Policies_k,h,nb_state, years,Q)

#@save "Vals_converter_120725.jld2" Vals_all
#@load "Vals_converter_120725.jld2" 
