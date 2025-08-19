#Vals_all, Policies_m, Policies_k = Bellman(1,11,1,h,12,P,1)

for s in 1:nb_state
    for t in 1:6
        @assert sum(Policies_k[s,1,t,:]) == 0 
    end
end

for q in 1:12
    for t in 1:6
        quota = q-1
        @assert sum(Policies_m[12,q,t,:]) >= 1.0 "Problem with (s = 12, q=$quota, T=$t)"
    end
end

