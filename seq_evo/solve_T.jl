#= using Pkg
Pkg.add("GLPK")
Pkg.add("JuMP")
=# 

#read in the father (X) and son (Y) matrices

using DelimitedFiles

X = readdlm("../results/seq_evo/sf_mat.txt")
Y = readdlm("../results/seq_evo/sons_mat.txt")

using JuMP, Ipopt , LinearAlgebra


#define model
m = Model(Ipopt.Optimizer)

# Set solver options (tolerance)
set_optimizer_attribute(m, "tol", 1e-15)
set_optimizer_attribute(m, "constr_viol_tol", 1e-15)

K = size(X,1)
N = size(X,2)

#define variables
@variable(m, 0<= T[1:K, 1:K] <= 1)
#@NLobjective(m, Min, sum(Y-T*X)^2 ), what I meant but Julia doesn't like it so code below.
@objective(m, Min, sum((Y[i, j] - sum(T[i, k] * X[k, j] for k in 1:K))^2 for i in 1:K, j in 1:N))
@constraint(m, [i = 1:K], sum(T[i, j] for j in 1:K) == 1)

#set starting value
# Set starting guess for T as a diagonal matrix of ones
for i in 1:K
    for j in 1:K
        if i == j
            set_start_value(T[i, j], 1.0)
        else
            set_start_value(T[i, j], 0.0)
        end
    end
end

# Solve the optimization problem
optimize!(m)

# Retrieve the optimal solution
optimal_T = value.(T)

#check columns sum to 1
sums = zeros(K)
for i in 1:K
    sums[i] = sum(optimal_T[:,i])
end

sums

#sum(optimal_T[:,2])

println("Optimal Solution for T:")
println(optimal_T)
