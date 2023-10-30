struct KleinGordon{T}
    d :: Int
    β₁ :: T
    β₂ :: T
end

#

function display_equation(eq::KleinGordon)
    print("    - β₁ = "); showfull(eq.β₁); println()
    print("    - β₂ = "); showfull(eq.β₂); println()
end

#

function N(eq::KleinGordon, v)
    v₁ = component(v, 1)

    N_ = zeros(eltype(v), Taylor(3order(space(space(v))))^1)

    project!(component(N_, 1), -v₁ + (eq.β₁ + eq.β₂*v₁)*v₁^2)

    return N_
end

function DN(eq::KleinGordon, v, n_T)
    v₁ = component(v, 1)

    DN_ = zeros(eltype(v), Taylor(n_T)^1, Taylor(n_T)^1)

    project!(component(DN_, 1, 1), Multiplication(-1 + (2eq.β₁ + 3eq.β₂*v₁)*v₁))

    return DN_
end

function opnorm_DN(eq::KleinGordon, v, X)
    v₁ = component(v, 1)
    return norm(-1 + (2eq.β₁ + 3eq.β₂*v₁)*v₁, X)
end

opnorm_D²N_abs(eq::KleinGordon, v, X, ϱ) = 2abs(eq.β₁) + 6abs(eq.β₂)*(norm(component(v, 1), X) + ϱ)

#

function f(eq::KleinGordon, w)
    w₁, w₂, w₃ = eachcomponent(w)

    f_ = zeros(eltype(w), Chebyshev(3order(space(space(w))))^3)

    project!(component(f_, 1), -w₁^2)
    project!(component(f_, 2), w₃)
    project!(component(f_, 3), -(eq.d-1)*w₁*w₃ + w₂ - (eq.β₁ + eq.β₂*w₂)*w₂^2)

    return f_
end

function Df(eq::KleinGordon, w, n_C)
    w₁, w₂, w₃ = eachcomponent(w)

    Df_ = zeros(eltype(w), Chebyshev(2order(space(space(w))) + n_C + 1)^3, Chebyshev(n_C + 1)^3)

    project!(component(Df_, 1, 1), Multiplication(-2w₁))

    project!(component(Df_, 2, 3), I)

    project!(component(Df_, 3, 1), Multiplication(-(eq.d-1)*w₃))
    project!(component(Df_, 3, 2), Multiplication(1 - (2eq.β₁ + 3eq.β₂*w₂)*w₂))
    project!(component(Df_, 3, 3), Multiplication(-(eq.d-1)*w₁))

    return Df_
end

function opnorm_Df(eq::KleinGordon, w, X)
    w₁, w₂, w₃ = eachcomponent(w)
    return max(2norm(w₁, X), 1, (eq.d-1)*(norm(w₁, X) + norm(w₃, X)) + norm(1 - (2eq.β₁ + 3eq.β₂*w₂)*w₂, X))
end

opnorm_D²f_abs(eq::KleinGordon, w, X, ϱ) = max(2, 2(eq.d-1) + 2abs(eq.β₁) + 6abs(eq.β₂)*(norm(component(w, 2), X) + ϱ))

#

ψ̂(eq::KleinGordon, c, Λ, Γ, μ, ℒ_y) =
    (abs(eq.β₁) + 3abs(eq.β₂) * (2c[1] + (1 + ℒ_y) * μ) / 2) * (1 + ℒ_y) * μ
