# ||=================================||
# || Stochastic Simulation Algorithm ||
# || Matthew Shin (Jan 2019)         ||
# ||=================================||

using Random
using Plots

# Define find_μ function
function find_μ(a, a0, r2)
    μ = 1
    sum = 0
    while r2*a0 > sum + a[μ]
        sum = sum + a[μ]
        μ += 1
    end
    return μ
end

# Define >> react function
function react(X, μ)
    if μ == 1
        X = X + [0, 1]
    elseif μ == 2
        X = X + [0, -2]
    end
    return X
end


# Step 0: Initialization
Nreact = 30000
c = [0.5, 0.005]  # >> reaction constants
X = [10, 3000]  # >> molecular population
Xn = X
t = 0.0  # time variable
tn = [t]
n = 0  # reaction counter
Random.seed!(0)

while n < Nreact
    global X, t, n

    # Step 1: Calculate reaction probabilities
    h = [X[1]*X[2], Int(X[2]*(X[2]-1)/2)]  # >> combinatorial function
    a = h .* c
    a0 = sum(a)

    # Step 2: Generate PRNs and calculate τ and μ
    r1 = rand()
    r2 = rand()
    τ = (1/a0)*log(1/r1)
    μ = find_μ(a, a0, r2)

    # Step 3:
    t = t + τ
    push!(tn, t)
    X = react(X, μ)
    append!(Xn, X)
    n += 1
end

# Plot simulation results
Xn = collect(reshape(Xn, (length(X), Nreact+1))')  # Reshape Xn to (Nreact+1) x N
tn = collect(tn)
plot(tn, Xn)
