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
        X = X + [0, 1, 0]
    elseif μ == 2
        X = X + [0, -1, 1]
    elseif μ == 3
        X = X + [0, 0, -2]
    end
    return X
end

# Step 0: Initialization
Nreact = 500000
prtbar = Int(Nreact/50)
c = [1.0, 0.01, 10.0]  # >> Reaction constants
X = [10, 1000, 1000]  # >> Molecular population
Xn = X
t = 0.0  # Time variable
tn = [t]
n = 0  # Reaction counter
Random.seed!(3)

print("\nCalculating...\n0% |")
while n < Nreact
    global X, t, n

    # Step 1: Calculate reaction probabilities
    h = [X[1]*X[2], X[2]*X[3], X[3]]  # >> Combinatorial function
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
    if n%prtbar == 0
        print('=')
    end
end
print("| 100%\n\n")

# Plot simulation results
Xn = collect(reshape(Xn, (length(X), Nreact+1))')  # Reshape Xn
Xn100 = Xn[1:100:end, :]  # Sample every 100 reactions
tn = collect(tn)
tn100 = tn[1:100:end]  # Sample every 100 reactions
println("Plotting...")
plot(tn100, Xn100)
#plot(Xn100[:, 2], Xn100[:, 3])
