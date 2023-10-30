function prove(filename::String)
    data = load(filename)

    eq, x, n_T_pad, n_C_pad, c, Λ, Γ, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ = data["equation"],
        data["x"], data["n_T_pad"], data["n_C_pad"], data["c"], data["Λ"], data["Γ"],
        data["l"], data["r_star"], data["L"], data["ℒ_x"], data["ℒ_y"], data["ν"], data["ϱ"]

    println("Starting the proof for the localized radial solution ", filename, " on ℝ^", eq.d, ":")
    display_data(eq, x, n_T_pad, n_C_pad, c, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ)

    return _prove(eq, x, n_T_pad, n_C_pad, c, Λ, Γ, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ)
end

function display_data(eq, x, n_T_pad, n_C_pad, c, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ)
    display_equation(eq)

    η = component(x, 1)
    v, w = component(x, 3), component(x, 4)

    n_T_num = order(space(space(v)))
    n_C_num = order(space(space(w)))

    δ = inv(l*Interval(r_star) + L)
    μ = norm(Interval.(η), Inf) + ϱ

    println("• Constant equilibrium:")
    print("    c = ["); foreach(cᵢ -> (showfull(cᵢ); print(", ")), c[1:end-1]); showfull(c[end]); println("]")

    println("• Approximate initial condition:")
    println("    ϕ̄ = $(coefficients(component(x, 2)))")

    println("• Approximate stable coordinate:")
    println("    η̄ = $(coefficients(η))")

    println("• Domain [-l, l] of the Taylor series: l = $l")
    println("  Numerical truncation order of the Taylor series: n_T_num = $n_T_num")
    println("  Domain [l r_star, l r_star + L] of the Chebyshev series: r_star = $r_star and L = $L")
    println("  Numerical truncation order of the Chebyshev series: n_C_num = $n_C_num")

    print("• Domain [0, δ] of the local graph of the centre-stable manifold in the centre direction: δ = "); showfull(δ); println()
    println("  Maximal error bound: ϱ = $ϱ")
    print("  Domain [-μ, μ] of the local graph of the centre-stable manifold in the stable direction: μ = "); showfull(μ); println()
    println("  Lipschitz constant of the local graph of the centre-stable manifold in the centre direction: ℒ_x = $ℒ_x")
    println("  Lipschitz constant of the local graph of the centre-stable manifold in the stable direction: ℒ_y = $ℒ_y")

    println("• Decay rate of the Chebyshev sequence space: ν = $ν")

    println("• Truncation order of the Taylor series for the computer-assisted proof: n_T = n_T_num + n_T_pad = $n_T_num + $n_T_pad")
    println("  Truncation order of the Chebyshev series for the computer-assisted proof: n_C = n_C_num + n_C_pad = $n_C_num + $n_C_pad")
end

_prove(eq::Union{KleinGordon{Interval{Float64}},FitzHughNagumo{Interval{Float64}}},
    x::Sequence{<:CartesianProduct,Vector{Float64}}, n_T_pad::Int, n_C_pad::Int,
    c::Sequence{CartesianPower{ParameterSpace},Vector{Interval{Float64}}},
    Λ::LinearOperator{CartesianPower{ParameterSpace},CartesianPower{ParameterSpace},Diagonal{Interval{Float64},Vector{Interval{Float64}}}},
    Γ::LinearOperator{CartesianPower{ParameterSpace},CartesianPower{ParameterSpace},Matrix{Interval{Float64}}},
    l::Float64, r_star::Float64, L::Float64, ℒ_x::Float64, ℒ_y::Float64, ν::Float64, ϱ::Float64) =
        __prove(eq, Interval.(x), n_T_pad, n_C_pad, c, Λ, Γ, Interval(l), Interval(r_star), Interval(L), Interval(ℒ_x), Interval(ℒ_y), Interval(ν), Interval(ϱ))

_prove(eq::SwiftHohenberg{Interval{Float64}},
    x::Sequence{<:CartesianProduct,Vector{Complex{Float64}}}, n_T_pad::Int, n_C_pad::Int,
    c::Sequence{CartesianPower{ParameterSpace},Vector{Interval{Float64}}},
    Λ::LinearOperator{CartesianPower{ParameterSpace},CartesianPower{ParameterSpace},Diagonal{Complex{Interval{Float64}},Vector{Complex{Interval{Float64}}}}},
    Γ::LinearOperator{CartesianPower{ParameterSpace},CartesianPower{ParameterSpace},Matrix{Complex{Interval{Float64}}}},
    l::Float64, r_star::Float64, L::Float64, ℒ_x::Float64, ℒ_y::Float64, ν::Float64, ϱ::Float64) =
        __prove(eq, Interval.(x), n_T_pad, n_C_pad, c, Λ, Γ, Interval(l), Interval(r_star), Interval(L), Interval(ℒ_x), Interval(ℒ_y), Interval(ν), Interval(ϱ))

function __prove(eq, x, n_T_pad, n_C_pad, c, Λ, Γ, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ)
    λ̂ = minimum(real(diag(coefficients(Λ))))

    η = component(x, 1)
    v, w = component(x, 3), component(x, 4)

    n_T_num = order(space(space(v)))
    n_T = n_T_num + n_T_pad

    n_C_num = order(space(space(w)))
    n_C = n_C_num + n_C_pad

    δ = inv(l*r_star + L)
    norm_η = norm(η, Inf)
    μ = norm_η + ϱ

    if isLPcontracting(eq, c, Λ, Γ, λ̂, δ, μ, ℒ_x, ℒ_y)
        println("-> Success: the Lyapunov-Perron operator is a contraction (cf. Proposition 2.2)")
        println("Computing Y, Z₁ and Z₂:")

        q = nspaces(space(v))

        #

        F_ = F(eq, x, c, Λ, Γ, l, r_star, L)

        DF_ = DF(eq, x, n_T, n_C, Λ, Γ, l, r_star, L)

        tail_N = copy(component(F_, 3))
        tail_f = copy(component(F_, 4))
        for i ∈ 1:q
            component(tail_N, i)[0:min(n_T-2, order(component(tail_N, i)))] .= 0
            component(tail_f, i)[0:min(n_C-1, order(component(tail_f, i)))] .= 0
            component(tail_f, i+q+1)[0:min(n_C-1, order(component(tail_f, i+q+1)))] .= 0
        end

        A = inv(mid.(project(DF_, codomain(DF_), codomain(DF_))))

        #

        X_T = ℓ¹()
        X_T_cartesian = NormedCartesianSpace(X_T, ℓ∞())

        X_C = ℓ¹(GeometricWeight(ν))
        X_C_cartesian = NormedCartesianSpace(X_C, ℓ∞())

        X = NormedCartesianSpace((ℓ∞(), ℓ∞(), X_T_cartesian, X_C_cartesian), ℓ∞())

        #

        Y = norm(A * F_, X) + max(
            l^2 * norm(tail_N, X_T_cartesian) / Interval((n_T + 1) * (n_T + eq.d - 1)),
            L * (ν + inv(ν)) * norm(tail_f, X_C_cartesian) / (4(n_C + 1))
            )
        print("    - Y = "); showfull(Y); println()

        n_C_max = Int((order(space(space(F_)[4]))-1)/n_C_num)*n_C + 1
        opnorm_A = opnorm(Interval.(A), X)
        opnorm_Df_ = opnorm_Df(eq, w, X_C)
        Z₁ = opnorm(A * DF_ - I, X) +
                max(
                    l^2 * opnorm_DN(eq, v, X_T) / Interval((n_T + 1) * (n_T + eq.d - 1)),
                    L * (ν + inv(ν)) * opnorm_Df_ / (4(n_C + 1))
                ) + opnorm_A *
                max(
                    2 / ν^(n_C_max + 1), # NOTE: 1 / ν^(n_C_max + 1) is also valid
                    r_star^(n_T + 1) * max(1, (n_T + 1) / l) + L * opnorm_Df_ / (2 * ν^(n_C + 2) * ((n_C + 2)^2 - 1))
                )
        print("    - Z₁ = "); showfull(Z₁); println()

        Z₂ = (opnorm_A + 1) *
            max(
                l^2 * opnorm_D²N_abs(eq, v, X_T, ϱ),
                L * (1 + ν) * opnorm_D²f_abs(eq, w, X_C, ϱ) / 2
                )
        print("    - Z₂ = "); showfull(Z₂); println()

        Q = LinearOperator(ParameterSpace()^q, ParameterSpace()^q × ParameterSpace()^q, [coefficients(Γ) ; coefficients(Γ*Λ)])
        opnorm_AQ = opnorm(component(A, :, 1:2) * Q, ℓ∞(), X)

        ρ = interval_of_existence(Y + opnorm_AQ*ℒ_y*norm_η, Z₁ + opnorm_AQ*ℒ_y, Z₂, ϱ)

        if !isempty(ρ)
            error = max(inf(ρ), inf(opnorm(Γ, Inf)*(ρ + (1 + ℒ_y)*μ)))
            println("-> Success: there exists a localized radial solution within a distance $error (in the C⁰-norm) of the numerical approximation (cf. Corollary 3.3)")
        else
            println("Failure: could not prove that there exists a localized radial solution nearby the numerical approximation (cf. Theorem 3.2)")
        end
    else
        println("Failure: could not prove that the Lyapunov-Perron operator is a contraction (cf. Proposition 2.2)")
    end
end

#

function F(eq, x, c, Λ, Γ, l, r_star, L)
    η, ϕ = component(x, 1), component(x, 2)
    v, w = component(x, 3), component(x, 4)
    q = nspaces(space(v))

    T = eltype(x)

    N_ = N(eq, v)
    f_ = f(eq, w)

    space_F₃ = Taylor(order(space(space(N_))) + 2)^q
    space_F₄ = Chebyshev(order(space(space(f_))) + 1)^(1 + 2q)
    space_F = ParameterSpace()^q × ParameterSpace()^q × space_F₃ × space_F₄
    F_ = zeros(T, space_F)

    n_T = order(space(space(v)))
    n_T_max = order(space(space_F₃))

    #

    component(F_, 1) .= component(w, 2:(1+q))(one(T)) .- c .- Γ * η
    component(F_, 2) .= component(w, (2+q):(1+2q))(one(T)) .+ Γ * Λ * η

    #

    F₃ = component(F_, 3)
    l² = pow(l, 2)
    for i ∈ 1:q
        component(F₃, i)[0] = component(v, i)[0] - ϕ[i]
        component(F₃, i)[1] = component(v, i)[1]
        for n ∈ 2:n_T_max
            component(F₃, i)[n] = l² * component(N_, i)[n-2]
            if n ≤ n_T
                component(F₃, i)[n] += n * (n + eq.d - 2) * component(v, i)[n]
            end
        end
    end

    #

    F₄ = component(F_, 4)
    integrate!(F₄, f_) .*= L/2
    rsub!(F₄, w)
    component(F₄, 1)[0] += inv(l*r_star)
    for i ∈ 1:q, n ∈ 0:n_T
        component(F₄, i+1)[0] += component(v, i)[n] * pow(r_star, n)
        component(F₄, i+q+1)[0] += n * component(v, i)[n] * pow(r_star, n) / l
    end

    #

    return F_
end

function DF(eq, x, n_T, n_C, Λ, Γ, l, r_star, L)
    v, w = component(x, 3), component(x, 4)
    q = nspaces(space(v))

    T = eltype(x)

    DN_ = DN(eq, v, n_T)
    Df_ = Df(eq, w, n_C)

    domain_DF = ParameterSpace()^q × ParameterSpace()^q × Taylor(n_T)^q × domain(Df_)
    space_x_padded = ParameterSpace()^q × ParameterSpace()^q × Taylor(n_T)^q × Chebyshev(n_C)^(1+2q)
    DF_ = zeros(T, domain_DF, space_x_padded)

    #

    component(DF_, 1, 1) .= (-).(Γ)
    mul!(component(DF_, 2, 1), Γ, Λ)
    D₄F₁ = component(DF_, 1, 4)
    D₄F₂ = component(DF_, 2, 4)
    for i ∈ 1:q
        component(D₄F₁, i, i+1) .= 2
        component(D₄F₁, i, i+1)[1,0] = 1
        component(D₄F₂, i, i+q+1) .= 2
        component(D₄F₂, i, i+q+1)[1,0] = 1
    end

    #

    D₂F₃ = component(DF_, 3, 2)
    D₃F₃ = component(DF_, 3, 3)
    l² = pow(l, 2)
    for j ∈ 1:q, i ∈ 1:q
        if i == j
            component(D₂F₃, i, i)[0,1] = -1
            component(D₃F₃, i, i)[0,0] = component(D₃F₃, i, i)[1,1] = 1
            for m ∈ 0:n_T, n ∈ 2:n_T
                component(D₃F₃, i, i)[n,m] = l² * component(DN_, i, i)[n-2,m]
                if n == m
                    component(D₃F₃, i, i)[n,n] += n * (n + eq.d - 2)
                end
            end
        else
            for m ∈ 0:n_T, n ∈ 2:n_T
                component(D₃F₃, i, j)[n,m] = l² * component(DN_, i, j)[n-2,m]
            end
        end
    end

    #

    rsub!(mul!(component(DF_, 4, 4), Integral(1), Df_, L/2, false), I)
    D₃F₄ = component(DF_, 4, 3)
    for i ∈ 1:q, n ∈ 0:n_T
        component(D₃F₄, i+1, i)[0,n] = pow(r_star, n)
        component(D₃F₄, i+q+1, i)[0,n] = n * pow(r_star, n) / l
    end

    #

    return DF_
end

#

function isLPcontracting(eq, c, Λ, Γ, λ̂, δ, μ, ℒ_x, ℒ_y)
    ψ̂_ = ψ̂(eq, c, Λ, Γ, μ, ℒ_y)

    condition_1 = sup(((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)) < inf(λ̂)

    condition_2 = inf(ℒ_y * (2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y))) ≥ sup(((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y))

    condition_3 = inf(ℒ_x) ≥ sup(
        (inv(2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)) +
        ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y) / ((2λ̂ - (3(eq.d - 1) * δ / 2 + 2ψ̂_) * (1 + ℒ_y)) * (2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)))) *
        ((eq.d - 1) * (1 + ℒ_y) / 2 + ((eq.d - 1) * δ / 2 + ψ̂_) * ℒ_x)
        )

    return condition_1 & condition_2 & condition_3
end
