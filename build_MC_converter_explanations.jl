MTBF = 588

p = (11/MTBF)/365 #gives a rough estimatin of p_i i+1

function proba(p,i,j)
    # j > i
    if (j-i)>5
        return 0
    else
        return binomial(1680-i+1,j-i)*(p^(j-i))*((1-p)^(1680-j+1))
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
    T = [0.0 for i in 1:12]

    for i in 11:-1:1
        T[i] = (1/(1 - proba_ii(p,i)))*(1 + sum(proba(p,i,j)*T[j] for j in (i+1):12))
    end

    return T[1]

end


nb_state = 12
nb_maint = 2



P_evac = vcat([100.0 for i in 1:11], [0.0]) # vector of capacity levels

# using dichotomy we see that we get the MTBF for p = 1.1341906625663963e-5
# binomial(1677, 8) overflows so we can't compute proba_i,i+n with n = 8. This is ok because this probability is very small, so we can neglect it. 
# proba(p,1,3) is 100 times smaller than proba(1,2)
# proba(p,1,2) is very close to proba(11,12)

#CCl: we take p(i,i+1) = p for each i, p(i,j) = 0 if j > i+1

P = zeros(nb_state, nb_state, nb_maint)

for i in 1:nb_state
    P[1,i,1] = 1.0 
end
    
for i in 1:(nb_state - 1)
    P[i,i,2] = 1.0 - p 
    P[i+1,i,2] = p
end
    
P[nb_state,nb_state,2] = 1.0

d = [11,0] # length of maintenance does not depend on state 

#100*(proba(p,1,2) - proba(p,11,12))/proba(p,11,12)

#tests

@assert all(x -> x <= 1, P)

tol = 1e-6  
@assert all(abs(sum(P[i, j, maint] for i in 1:nb_state) - 1) ≤ tol
            for j in 1:nb_state, maint in 1:2)