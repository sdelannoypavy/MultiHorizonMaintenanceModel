using Gurobi
using JuMP  
using Distributions
using Random
using StatsBase


include("parameters_weakly_coupled.jl")


"""
---------------------------------------------------------------------------------------------
Construct model
---------------------------------------------------------------------------------------------
"""

function create_α(x0,k0)

    α = [0.0 for c in 1:C, x in 1:n_max, k in 0:Q]

    for c in 1:C
        α[c,x0[c],k0+1] = 1
    end

    return α

end

α = create_α(x0,k0)

println("started optimization")
model = Model(Gurobi.Optimizer)
#set_optimizer_attribute(model, "OutputFlag", 0)


@variable(model, q[c=1:C, x=1:n[c], k in 0:Q, m in 0:1, κ in 0:Q, t in 1:T] >= 0)     
@variable(model, A[κ in 0:Q, t in 1:T] >= 0)   
@variable(model, μ[k in 0:Q, t in 1:T] >= 0)  

M = 100.0

@objective(model, Min, M*sum(δ[c,x,m+1,κ+1,t]*q[c,x,k,m,κ,t] for c in 1:C, x in 1:n[c], k in 0:Q, m in 0:1, κ in 0:Q, t in 1:T))

# do we really want one vector of probability transitions for each t?
@constraint(model, [c in 1:C, t in 2:T, x in 1:n[c], k in 0:Q], sum(q[c,x,k,m,κ,t] for m in 0:1, κ in 0:Q) 
    == sum(p[c,x,k+1,x′,k′+1,m+1,κ+1,t]*q[c,x,k,m,κ,t-1] for x′ in 1:n[c], k′ in 0:Q, m in 0:1, κ in 0:Q))

@constraint(model, [c in 1:C, m in 0:1, κ in 0:Q, t in 1:T], sum(q[c,x,k,m,κ,t] for x in 1:n[c], k in 0:Q, m in 0:1) == A[κ,t])

@constraint(model, [c in 1:C, x in 1:n[c], k in 0:Q], sum(q[c,x,k,m,κ,1] for m in 0:1, κ in 0:Q) == α[c,x,k+1])

@constraint(model, [c in 1:C, x in 1:n[c], k in 0:(Q-1), m in 0:1, κ in (k+1):Q, t in 1:T], q[c,x,k,m,κ,t] == 0)

@constraint(model, [c in 1:C, k in 0:Q, t in 1:T], sum(q[c,x,k,m,κ,t] for x in 1:n[c], m in 0:1, κ in 0:Q) == μ[k,t])


"""
---------------------------------------------------------------------------------------------
Optimize
---------------------------------------------------------------------------------------------
"""

optimize!(model)

status = termination_status(model)
cost = objective_value(model)

A = [value(A[κ,1]) for κ in 0:Q]
A_opt = argmax(A) - 1
m_opt = [argmax([sum(value(q[c,x,k,m,κ,1]) for x in 1:n[c], k in 0:Q, κ in 0:Q) for m in 0:1]) - 1 for c in 1:C]


if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
    println("total cost: ",cost)
    println("optimal κ at step 1: ", A_opt)
    println("optimal m at step 1: ",  m_opt)
else
    println("Aucune solution optimale trouvée.")
    println(status)
end    


