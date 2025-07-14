MTBF = 588

psim = (12/MTBF)/365 #gives a rough estimatin of p_i i+1

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

P = vcat([100.0 for i in 1:11], [0.0]) # vector of capacity levels

# using dichotomy we see that we get the MTBF for p = 1.1341906625663963e-5
# binomial(1677, 8) overflows so we can't compute proba_i,i+n with n = 8. This is ok because this probability is very small, so we can neglect it. 
# proba(p,1,3) is 100 times smaller than proba(1,2)
# proba(p,1,2) is very close to proba(11,12)

#CCl: we take p(i,i+1) = psim for each i, p(i,j) = 0 if j > i+1