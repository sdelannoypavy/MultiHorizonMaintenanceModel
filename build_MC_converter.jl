true_MTBF_converter = 588 # in years

nb_state = 12
P_evac = vcat([100.0 for i in 1:11], [0.0]) # vector of capacity levels

function build_P_converter(MTBF,nb_state)

    p = (11/MTBF)/365

    P = zeros(nb_state, nb_state, 2)

    for i in 1:nb_state
        P[1,i,1] = 1.0 
    end
        
    for i in 1:(nb_state - 1)
        P[i,i,2] = 1.0 - p 
        P[i+1,i,2] = p
    end
        
    P[nb_state,nb_state,2] = 1.0

    return P
end

#P = build_P_converter(true_MTBF)