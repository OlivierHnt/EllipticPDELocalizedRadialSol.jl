struct FitzHughNagumo{T}
    d :: Int
    ϵ :: T
    β₁ :: T
    β₂ :: T
    β₃ :: T
    β₄ :: T
end

#

function display_equation(eq::FitzHughNagumo)
    print("    - ϵ = "); showfull(eq.ϵ); println()
    print("    - β₁ = "); showfull(eq.β₁); println()
    print("    - β₂ = "); showfull(eq.β₂); println()
    print("    - β₃ = "); showfull(eq.β₃); println()
    print("    - β₄ = "); showfull(eq.β₄); println()
end

#

function N(eq::FitzHughNagumo, v)
    v₁, v₂, v₃ = eachcomponent(v)

    space_N = Taylor(3order(space(space(v))))^3
    N_ = Sequence(space_N, zeros(eltype(v), dimension(space_N)))

    ϵ⁻¹ = inv(eq.ϵ)

    project!(component(N_, 1), ϵ⁻¹^2 * (1 - v₁^2)*v₁ - ϵ⁻¹ * (eq.β₂*v₂ + eq.β₃*v₃ + eq.β₁))
    project!(component(N_, 2), v₁ - v₂)
    project!(component(N_, 3), (v₁ - v₃)/eq.β₄^2)

    return N_
end

function DN(eq::FitzHughNagumo, v)
    v₁, v₂, v₃ = eachcomponent(v)

    DN_ = LinearOperator(space(v), space(v), zeros(eltype(v), length(v), length(v)))

    ϵ⁻¹ = inv(eq.ϵ)
    β₄⁻² = inv(eq.β₄^2)

    project!(component(DN_, 1, 1), Multiplication(ϵ⁻¹^2 * (1 - 3v₁^2)))
    project!(component(DN_, 1, 2), -ϵ⁻¹ * eq.β₂ * I)
    project!(component(DN_, 1, 3), -ϵ⁻¹ * eq.β₃ * I)


    project!(component(DN_, 2, 1), I)
    project!(component(DN_, 2, 2), -I)


    project!(component(DN_, 3, 1), β₄⁻² * I)
    project!(component(DN_, 3, 3), -β₄⁻² * I)

    return DN_
end

function opnorm_DN(eq::FitzHughNagumo, v, X)
    v₁ = component(v, 1)

    ϵ⁻¹ = inv(eq.ϵ)

    return max(norm(ϵ⁻¹^2 * (1 - 3v₁^2), X) + abs(ϵ⁻¹) * (abs(eq.β₂) + abs(eq.β₃)), 2, 2/eq.β₄^2)
end

opnorm_D²N_abs(eq::FitzHughNagumo, v, X, ϱ) = 6inv(abs(eq.ϵ^2))*(norm(component(v, 1), X) + ϱ)

#

function f(eq::FitzHughNagumo, w)
    w₁, w₂, w₃, w₄, w₅, w₆, w₇ = eachcomponent(w)

    space_f = Chebyshev(3order(space(space(w))))^7
    f_ = Sequence(space_f, zeros(eltype(w), dimension(space_f)))

    ϵ⁻¹ = inv(eq.ϵ)

    project!(component(f_, 1), -w₁^2)
    project!(component(f_, 2), w₅)
    project!(component(f_, 3), w₆)
    project!(component(f_, 4), w₇)
    project!(component(f_, 5), -(eq.d-1)*w₁*w₅ - (ϵ⁻¹^2 * (1 - w₂^2)*w₂ - ϵ⁻¹ * (eq.β₂*w₃ + eq.β₃*w₄ + eq.β₁)))
    project!(component(f_, 6), -(eq.d-1)*w₁*w₆ - (w₂ - w₃))
    project!(component(f_, 7), -(eq.d-1)*w₁*w₇ - (w₂ - w₄)/eq.β₄^2)

    return f_
end

function Df(eq::FitzHughNagumo, w)
    w₁, w₂, w₃, w₄, w₅, w₆, w₇ = eachcomponent(w)

    n_C = order(space(space(w)))
    domain_Df = Chebyshev(3n_C + 1)^7
    codomain_Df = Chebyshev(n_C + 1)^7
    Df_ = LinearOperator(domain_Df, codomain_Df, zeros(eltype(w), dimension(codomain_Df), dimension(domain_Df)))

    ϵ⁻¹ = inv(eq.ϵ)
    β₄⁻² = inv(eq.β₄^2)

    project!(component(Df_, 1, 1), Multiplication(-2w₁))

    project!(component(Df_, 2, 5), I)

    project!(component(Df_, 3, 6), I)

    project!(component(Df_, 4, 7), I)

    project!(component(Df_, 5, 1), Multiplication(-(eq.d-1)*w₅))
    project!(component(Df_, 5, 2), Multiplication(-ϵ⁻¹^2 * (1 - 3w₂^2)))
    project!(component(Df_, 5, 3), ϵ⁻¹ * eq.β₂ * I)
    project!(component(Df_, 5, 4), ϵ⁻¹ * eq.β₃ * I)
    project!(component(Df_, 5, 5), Multiplication(-(eq.d-1)*w₁))

    project!(component(Df_, 6, 1), Multiplication(-(eq.d-1)*w₆))
    project!(component(Df_, 6, 2), -I)
    project!(component(Df_, 6, 3), I)
    project!(component(Df_, 6, 6), Multiplication(-(eq.d-1)*w₁))

    project!(component(Df_, 7, 1), Multiplication(-(eq.d-1)*w₇))
    project!(component(Df_, 7, 2), -β₄⁻² * I)
    project!(component(Df_, 7, 4), β₄⁻² * I)
    project!(component(Df_, 7, 7), Multiplication(-(eq.d-1)*w₁))

    return Df_
end

function opnorm_Df(eq::FitzHughNagumo, w, X)
    w₁, w₂, w₃, w₄, w₅, w₆, w₇ = eachcomponent(w)

    ϵ⁻¹ = inv(eq.ϵ)

    return max(2norm(w₁, X), 1,
        (eq.d-1)*(norm(w₁, X) + norm(w₅, X)) + norm(ϵ⁻¹^2 * (1 - 3w₅^2), X) + abs(ϵ⁻¹) * (abs(eq.β₂) + abs(eq.β₃)),
        (eq.d-1)*(norm(w₁, X) + norm(w₆, X)) + 2,
        (eq.d-1)*(norm(w₁, X) + norm(w₇, X)) + 2/eq.β₄^2)
end

opnorm_D²f_abs(eq::FitzHughNagumo, w, X, ϱ) = max(2, 2(eq.d-1) + 6inv(abs(eq.ϵ^2))*(norm(component(w, 2), X) + ϱ))

#

ψ̂(eq::FitzHughNagumo, c, Λ, Γ, μ, ℒ_y) =
    (3inv(abs(eq.ϵ^2)) * (2c[1] + (1 + ℒ_y) * μ) * (1 + ℒ_y) * μ / 2) * opnorm(inv(Γ * Λ), Inf) * opnorm(Γ, Inf)
