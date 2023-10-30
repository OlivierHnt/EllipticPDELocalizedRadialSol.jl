struct SwiftHohenberg{T}
    d :: Int
    β₁ :: T
    β₂ :: T
    β₃ :: T
    β₄ :: T
end

#

function display_equation(eq::SwiftHohenberg)
    print("    - β₁ = "); showfull(eq.β₁); println()
    print("    - β₂ = "); showfull(eq.β₂); println()
    print("    - β₃ = "); showfull(eq.β₃); println()
    print("    - β₄ = "); showfull(eq.β₄); println()
end

#

function N(eq::SwiftHohenberg, v)
    v₁, v₂ = eachcomponent(v)

    N_ = zeros(eltype(v), Taylor(3order(space(space(v))))^2)

    project!(component(N_, 1), eq.β₄*v₁ - v₂)
    project!(component(N_, 2), eq.β₄*v₂ - eq.β₁*v₁ - (eq.β₂ + eq.β₃*v₁)*v₁^2)

    return N_
end

function DN(eq::SwiftHohenberg, v, n_T)
    v₁, v₂ = eachcomponent(v)

    DN_ = zeros(eltype(v), Taylor(n_T)^2, Taylor(n_T)^2)

    project!(component(DN_, 1, 1), eq.β₄ * I)
    project!(component(DN_, 1, 2), -I)

    project!(component(DN_, 2, 1), Multiplication(-eq.β₁ - (2eq.β₂ + 3eq.β₃*v₁)*v₁))
    project!(component(DN_, 2, 2), eq.β₄ * I)

    return DN_
end


function opnorm_DN(eq::SwiftHohenberg, v, X)
    v₁ = component(v, 1)
    return max(abs(eq.β₄) + 1, abs(eq.β₄) + norm(eq.β₁ + (2eq.β₂ + 3eq.β₃*v₁)*v₁, X))
end

opnorm_D²N_abs(eq::SwiftHohenberg, v, X, ϱ) = 2abs(eq.β₂) + 6abs(eq.β₃)*(norm(component(v, 1), X) + ϱ)

#

function f(eq::SwiftHohenberg, w)
    w₁, w₂, w₃, w₄, w₅ = eachcomponent(w)

    f_ = zeros(eltype(w), Chebyshev(3order(space(space(w))))^5)

    project!(component(f_, 1), -w₁^2)
    project!(component(f_, 2), w₄)
    project!(component(f_, 3), w₅)
    project!(component(f_, 4), -(eq.d-1)*w₁*w₄ - eq.β₄*w₂ + w₃)
    project!(component(f_, 5), -(eq.d-1)*w₁*w₅ - eq.β₄*w₃ + eq.β₁*w₂ + (eq.β₂ + eq.β₃*w₂)*w₂^2)

    return f_
end

function Df(eq::SwiftHohenberg, w, n_C)
    w₁, w₂, w₃, w₄, w₅ = eachcomponent(w)

    Df_ = zeros(eltype(w), Chebyshev(2order(space(space(w))) + n_C + 1)^5, Chebyshev(n_C + 1)^5)

    project!(component(Df_, 1, 1), Multiplication(-2w₁))

    project!(component(Df_, 2, 4), I)

    project!(component(Df_, 3, 5), I)

    project!(component(Df_, 4, 1), Multiplication(-(eq.d-1)*w₄))
    project!(component(Df_, 4, 2), -eq.β₄ * I)
    project!(component(Df_, 4, 3), I)
    project!(component(Df_, 4, 4), Multiplication(-(eq.d-1)*w₁))

    project!(component(Df_, 5, 1), Multiplication(-(eq.d-1)*w₅))
    project!(component(Df_, 5, 2), Multiplication(eq.β₁ + (2eq.β₂ + 3eq.β₃*w₂)*w₂))
    project!(component(Df_, 5, 3), -eq.β₄ * I)
    project!(component(Df_, 5, 5), Multiplication(-(eq.d-1)*w₁))

    return Df_
end

function opnorm_Df(eq::SwiftHohenberg, w, X)
    w₁, w₂, w₃, w₄, w₅ = eachcomponent(w)
    return max(2norm(w₁, X), 1,
        (eq.d-1)*(norm(w₁, X) + norm(w₄, X)) + abs(eq.β₄) + 1,
        (eq.d-1)*(norm(w₁, X) + norm(w₅, X)) + abs(eq.β₄) + norm(eq.β₁ + (2eq.β₂ + 3eq.β₃*w₂)*w₂, X))
end

opnorm_D²f_abs(eq::SwiftHohenberg, w, X, ϱ) = max(2, 2(eq.d-1) + 2abs(eq.β₂) + 6abs(eq.β₃)*(norm(component(w, 2), X) + ϱ))

#

ψ̂(eq::SwiftHohenberg, c, Λ, Γ, μ, ℒ_y) =
    (abs(eq.β₂) + 3abs(eq.β₃) * (2c[2] + (1 + ℒ_y) * μ) / 2) * (1 + ℒ_y) * μ * opnorm(inv(Γ * Λ), Inf) * opnorm(Γ, Inf)
