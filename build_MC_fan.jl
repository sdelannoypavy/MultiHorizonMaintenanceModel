true_MTBF = 40 # in years

nb_state = 8
P_evac = [100.0, 100.0, 100.0, 90.0, 80.0, 70.0, 60.0, 50.0] # vector of capacity levels

function build_P_fan(MTBF)

    p = (7/MTBF)/365

    P = zeros(nb_state, nb_state, 2)

    for i in 1:nb_state
        P[1,i,1] = 1.0 
    end
        
    for i in 1:(nb_state - 1)
        sum = 0
        for j in (i+1):nb_state
            P[j,i,2] = binomial(12-i+1,j-i)*(p^(j-i))*((1-p)^(12-j+1))
            sum += P[j,i,2]
        end
        P[i,i,2] = 1.0 - sum
    end
        
    P[nb_state,nb_state,2] = 1.0

    return P
end

P = build_P_fan(true_MTBF)

