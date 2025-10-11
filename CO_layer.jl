function CO_layer(θ; H, x, k, L, Q, nb_scen)

    println(x,k,H)

    θx = [θ[c] for c in 1:C]
    θk = θ[C+1]

    L = Int(round(L))

    c = zeros(L, Q+1, 30)

    T_per = 62 - count(==( -1 ), h[H,1,:]) 

    for m in 1:L
        for κ in 0:Q
            for t in 1:30
                total_cost = 0
                
                main_dec = [0 for i in 1:C]
                for c in 1:C
                    if d[c] <= length_m[m]
                        main_dec[c] = 1
                    end
                end

                Policies_m_t_list = [create_schedule_m(main_dec[c],t,T_per) for c in 1:C]
                Policies_κ_t = create_schedule_κ(κ,t,T_per)

                for scen in nb_scen
                    cost, _ = simulate_period(Policies_m_t_list,Policies_κ_t,x,h[H,scen,:],pr[H,:],n,d,P_daily)
                    total_cost += (1/nb_scen) * cost
                end
                c[m,κ+1,t] = total_cost
            end
        end
    end

    
    model = Model(Gurobi.Optimizer)
    @variable(model, δ[1:L, 0:k, 1:30], Bin)
    @variable(model, x_end[c=1:C, x=1:n[c]] >= 0)
    @variable(model, m[c=1:C], Bin)
    @variable(model, k_end, Bin)
    @variable(model, leng, Bin)

    @constraint(model, sum(δ[l, κ, t] for l in 1:L, κ in 0:k, t in 1:30) == 1)

    @objective(model, Min, sum(c[l,κ+1,t]*δ[l,κ,t] for l in 1:L, κ in 0:k, t in 1:30) + sum(θx[c]*sum(s*x_end[c,s] for s in 1:n[c]) for c in 1:C) + θk*k_end)

    @constraint(model,k_end == k - sum(κ*sum(δ[l, κ, t] for l in 1:L, t in 1:30) for κ in 0:k))
    @constraint(model,leng == sum(length_m[l]*δ[l, κ, t] for l in 1:L, κ in 0:k, t in 1:30))

    @constraint(model, [c in 1:C], m[c] <= 100*(leng - d[c] + 1.0))

    @constraint(model, [c in 1:C, s in 1:n[c]], !m[c] --> {x_end[c,s] == P_period[c][s,x[c]]})
    
    @constraint(model, [c in 1:C, s in 1:(n[c]-1)], m[c] --> {x_end[c,s] == 0.0})
    @constraint(model, [c in 1:C], m[c] --> {x_end[c,n[c]] == 1.0})


    try 
        optimize!(model)
        #return(cost, κ, m, t, l_opt)
        κ = sum(kk * value(δ[l, kk, t]) for l in 1:L, kk in 0:k, t in 1:30; init=0)
        #t = sum(t*value(δ[l, κ, t]) for l in 1:L, κ in 0:k, t in 1:30)
        #l_opt = sum(value(l))

        m = [value(m[c]) for c in 1:C]
        return([m; κ])
    catch e # error when no failure
        return([[0 for c in 1:C]; 0])
    end

end

