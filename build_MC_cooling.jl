MTBF_tot = 43

p = (12/MTBF)/365 

nb_state_water = 3 #could be more

MTBF_water = MTBF_tot # could be something else
p_water = (2/MTBF_water)/365 
P_water = [1 - p_water  0                 0
           p_water      1 - p_water       0
           0            p_water           1]

p_pump = 9.55e-5 # gives MTBF_pump around 43
P_pump = [1 - 2*p_pump*(1- p_pump)  0                 0
          2*p_pump*(1- p_pump)      1 - p_pump        0
          0                         p_pump            1]

MTBF_pump = ((1/(2*p_pump*(1- p_pump))) + 1/p_pump)/365

p_fan = 5.224e-5 #gives MTBF around 43
function proba(p,i,j)
    # j > i
    if (j-i)>5
        return 0
    else
        return binomial(12-i+1,j-i)*(p^(j-i))*((1-p)^(12-j+1))
    end
end 

function proba_ii(p,i)
    p_ii = 1
    for j in (i+1):12
        p_ii -= proba(p,i,j)
    end
    return p_ii
end

function compute_MTBF(p)
    # p is the probability of failure during one day, considering one single sub module

    # T = vector of time before absorption
    T = [0.0 for i in 1:8]

    for i in 7:-1:1
        T[i] = (1/(1 - proba_ii(p,i)))*(1 + sum(proba(p,i,j)*T[j] for j in (i+1):8))
    end

    return T[1]/365

end

P_fan = zeros(8,8)
for i in 1:7
    P_fan[i,i] = 1 - proba(p_fan,i,i+1)
    P_fan[i+1,i] = proba(p_fan,i,i+1)
end
P_fan[8,8] = 1.0

# We set to zero the probability of losing two fans simultaneously

nb_state = nb_state_water*3*8
nb_maint = 8 # 2³

P_no_maintenance = zeros(nb_state, nb_state)
P_m_water = zeros(nb_state, nb_state)
P_m_pump = zeros(nb_state, nb_state)
P_m_fan = zeros(nb_state, nb_state)

P_evac = [100 for s in 1:nb_state]

for state_water in 1:nb_state_water
    for state_pump in 1:3
        for state_fan in 1:8

            j = 24 * (state_water-1) + 8*(state_pump - 1) + state_fan

            P_m_water[8*(state_pump - 1) + state_fan,j] = 1.0
            P_m_pump[24 * (state_water-1) + state_fan,j] = 1.0
            P_m_fan[24 * (state_water-1) + 8*(state_pump - 1) + 1,j] = 1.0

            if (state_water == nb_state_water) || (state_pump == 3) || (state_fan == 8)
                P_evac[j] = 0
            end

            for state_water_target in 1:nb_state_water
                for state_pump_target in 1:3
                    for state_fan_target in 1:8
                        i = 24 * (state_water_target-1) + 8*(state_pump_target - 1) + state_fan_target
                        P_no_maintenance[i,j] = P_water[state_water_target,state_water]*P_pump[state_pump_target,state_pump]*P_fan[state_fan_target,state_fan]
                    end
                end
            end
        end
    end
end

d =  vcat([1 for i in 1:7], [0]) # we need two days to change the 2 pumps, but for now we assume that every maintenance takes 1 day


P = repeat(P_no_maintenance, 1, 1, nb_maint)

for m_water in 0:1 # 0 means do maintenance
    for m_pump in 0:1
        for m_fan in 0:1

            P_res = P_no_maintenance

            if m_water == 0
                P_res *= P_m_water
            end

            if m_pump == 0
                P_res *= P_m_pump
            end

            if m_fan == 0
                P_res *= P_m_fan
            end

            P[:,:,4 * (m_water) + 2*(m_pump) + m_fan + 1] = P_res

        end
    end
end


# tests


@assert all(x -> x <= 1, P)

tol = 1e-6  
@assert all(abs(sum(P[i, j, maint] for i in 1:nb_state) - 1) ≤ tol
            for j in 1:nb_state, maint in 1:nb_maint)

@assert all(abs(sum(P_water[i, j] for i in 1:nb_state_water) - 1) ≤ tol
        for j in 1:3)

@assert all(abs(sum(P_water[i, j] for i in 1:3) - 1) ≤ tol
        for j in 1:3)

@assert all(abs(sum(P_fan[i, j] for i in 1:8) - 1) ≤ tol
        for j in 1:3)

function test_sums()
    for i in axes(P, 2)
        for j in axes(P, 3)
            s = sum(P[n, i, j] for n in axes(P, 1))
            if abs(s - 1) > tol
                println("⚠️  Somme incorrecte pour i = $i, j = $j : somme = $s")
            end
        end
    end
end

