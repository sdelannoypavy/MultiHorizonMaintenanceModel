include("parameters_weakly_coupled.jl")
include("weakly_coupled3.jl")
include("benchmark_policies.jl")

using StatsPlots, CSV, StatsBase

cache = Dict{Tuple{Tuple{Vararg{Int64}}, Int64, Int64}, Tuple{Vector{Int64}, Int64}}()

function simulate_period(Policies_m_t_list,Policies_κ_t,last_state,h,pr,n,d,P_daily)

        # simulate states only over a strategic period for one failure scenario, weather scenario given by h, starting from state

        T_per = 62 - count(==( -1 ), h)

        u = [0 for c in 1:C, i in 1:T_per]
        q = [0 for c in 1:C, i in 1:T_per]
    
        for c in 1:C
            q[c,1] = d[c]*round(Int,Policies_m_t_list[c][1])
        end

        for c in 1:C
            if (q[c,1]>0)&&(h[1]==1)
                u[c,1] = 1
            end
        end
    
        costs = []
        ENE = []
        indispo = []
        new_cost = (1 - Policies_κ_t[1]) * max(pr[1] - Capacity(last_state,n,C,P_evac,q[:,1]), 0)
        new_ENE = max(pr[1] - Capacity(last_state,n,C,P_evac,q[:,1]), 0)
        if new_ENE > 0
            new_indispo = 1
        else
            new_indispo = 0
        end

        push!(costs, new_cost)
        push!(ENE, new_ENE)
        push!(indispo, new_indispo)
    
        for t in 1:T_per-1

            new_state = [1 for i in 1:C]
    
            for c in 1:C
                if u[c,t] == 0 #no maintenance
                    probas = P_daily[c][:, last_state[c]]  
                    new_state[c] = sample(1:n[c], Weights(probas))     
                else       
                    new_state[c] = 1
                end
            end
        
    
            for c in 1:C
                q[c,t+1] = q[c,t] - u[c,t] + d[c]*round(Int,Policies_m_t_list[c][t+1])
                if (q[c,t+1]>0)&&(h[t+1]==1)
                    u[c,t+1] = 1
                end
            end

            new_cost = (1 - Policies_κ_t[t+1]) * max(pr[t+1] - Capacity(new_state,n,C,P_evac,q[:,t+1]),0)
            new_ENE = max(pr[t+1] - Capacity(new_state,n,C,P_evac,q[:,t+1]),0)

            if new_ENE > 0
                new_indispo = 1
            else
                new_indispo = 0
            end
    
            last_state = new_state
            push!(costs, new_cost)
            push!(ENE, new_ENE)
            push!(indispo, new_indispo)
        end
    
        for c in 1:C
            if u[c,T_per] == 0 #no maintenance
                probas = P_daily[c][:, last_state[c]]   
                last_state[c] = sample(1:n[c], Weights(probas))     
            else       
                last_state[c] = 1
            end
        end

    return sum(costs), last_state, sum(ENE), sum(indispo)

end 


function simulate_concession_period(H,h_sim, pr, n, years, Q, d, P_daily,file_name_simu,file_name_policies,method,δ,p,display::Bool, C)

    # simulate stationary policies

    # simulate states over the entire concession period

    total_cost = 0
    total_ENE = 0
    total_indispo = 0

    last_state = [1 for c in 1:C]#at the beginning of the concession period, all components are new
    k = Q #we have not use any free maintenance day yet


    for t in 1:(6*years)

        if display
            println("Simulating cost for t = ",t," with state = ",last_state," and k = ", k)
        end

        if (method == "fluid") || (method == "fluid_no_w")
            Policies_m_t_list, Policies_κ_t = π(H, last_state, k, t, h_sim[t,:], file_name_policies, δ, p, pr, n, d, P_daily)
        elseif method == "fluid2"
            Policies_m_t_list, Policies_κ_t = π2(H, last_state, k, t, h_sim[t,:], file_name_policies, δ, p, pr, n, d, P_daily)
        elseif method == "benchmark1"
            Policies_m_t_list, Policies_κ_t = π_benchmark1(last_state,k,t)
        elseif method == "bellman"
            Policies_m_t_list, Policies_κ_t = π_bellman(last_state,k,t)
        else method == "learning"
            Policies_m_t_list, Policies_κ_t = π_learning(last_state,k,t,h_sim[t,:],file_name_policies,δ,p,pr,n,d,P_daily)
        end

        m_opti = [sum(Policies_m_t_list[c]) for c in 1:C]
        κ_opti = sum(Policies_κ_t) 

        new_cost, new_state, new_ENE, new_indispo = simulate_period(Policies_m_t_list,Policies_κ_t,last_state,h_sim[t,:],pr[t,:],n,d,P_daily)

        row = DataFrame(reshape([last_state; k; t; m_opti; κ_opti; new_cost; new_ENE; new_indispo], 1, :), :auto) 
        CSV.write(file_name_simu, row; append=true, writeheader=false) 

        if t%6 == 0
            k = Q #end of the year
        else
            k -= κ_opti
        end
        
        total_cost += new_cost
        total_ENE += new_ENE
        total_indispo += new_indispo
        last_state = new_state

    end

    return total_cost, total_ENE, total_indispo

end



function simulate(nb_w,nb_d,H,h,pr,n,years,Q,d,P_daily,file_name_simu,file_name_policies,method,δ,p,cost_file,display::Bool,C)

    for scen_w in 1:nb_w

        cost_mean = 0
        ENE_mean = 0
        indispo_mean = 0

        h_sim = [0 for T in 1:180, t in 1:62]

        for T in 1:180
            s = rand(1:80)
            h_sim[T,:] = h[T,s,:]
        end

        for _ in 1:nb_d
            new_cost, new_ENE, new_indispo = simulate_concession_period(H,h_sim, pr, n, years, Q, d, P_daily, file_name_simu, file_name_policies, method, δ, p, display, C)
            cost_mean += new_cost
            ENE_mean += new_ENE
            indispo_mean += new_indispo

            #row = DataFrame(reshape([method; new_cost], 1, :), :auto) 
        end

        cost_mean = cost_mean/nb_d
        ENE_mean = ENE_mean/nb_d
        indispo_mean = indispo_mean/nb_d

        row = DataFrame(reshape([method; cost_mean; ENE_mean; indispo_mean], 1, :), :auto) 
        CSV.write(cost_file, row; append=true, writeheader=false) 

    end


end

function collect_policies(file_name, C::Int64)

    df = CSV.read(file_name, DataFrame)

    for i in 1:nrow(df)
        row = df[i, :]
        vect_row = collect(row)
        
        x0 = vect_row[1:C]
        k0 = vect_row[C+1]
        t = vect_row[C+2]
        m_opt = vect_row[(C+3):(2C+2)]
        κ_opt = vect_row[2C+3]

        cle = (Tuple(x0), k0, t)
        cache[cle] = (m_opt, κ_opt)
        
    end

end
