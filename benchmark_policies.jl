# greedy policy

function π(x0,k0,t,h,file_name_policies)

    m_opt = [0 for c in 1:C]
    κ_opt = 0

    if t%6 == 4 #summer

        for c in 1:C
            if d[c] <= k0
                m_opt[c] = 1
            end 
        end

        κ_opt = k0
    end

    for c in 1:C
        if x0[c] == n[c] #failure 
            m_opt[c] = 1
            κ_opt = k0
        end
    end

    Policies_m_t_list = [create_schedule_m(m_opt[c],t_start,62) for c in 1:C]
    Policies_κ_t = create_schedule_κ(κ_opt,t_start,62)

    return Policies_m_t_list, Policies_κ_t

end