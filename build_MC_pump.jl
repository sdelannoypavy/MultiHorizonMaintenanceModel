nb_state = 3
P_evac = [100.0, 100.0, 0.0]

p_pump = 9.55e-5 # gives MTBF_pump around 43

P_pump = [1 - 2*p_pump*(1- p_pump)  0                 0
          2*p_pump*(1- p_pump)      1 - p_pump        0
          0                         p_pump            1]



function build_P_pump(p_pump)
        
    P = zeros(nb_state, nb_state, 2)
        
    for i in 1:nb_state
        P[1,i,1] = 1.0 
    end

    P[1,1,2] = 1 - 2*p_pump*(1- p_pump)
    P[2,1,2] = 2*p_pump*(1- p_pump)  
    P[2,2,2] = 1 - p_pump
    P[3,2,2] = p_pump 
    P[3,3,2] = 1          

    return P

end

true_MTBF = ((1/(2*p_pump*(1- p_pump))) + 1/p_pump)/365

P = build_P_pump(p_pump)