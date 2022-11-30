module EllipticPDELocalizedRadialSol

    import JLD2: load
    import LinearAlgebra: Diagonal, diag
    using RadiiPolynomial

include("klein_gordon.jl")
    export KleinGordon

include("swift_hohenberg.jl")
    export SwiftHohenberg

include("fitzhugh_nagumo.jl")
    export FitzHughNagumo

include("proof.jl")
    export prove

function julia_main()::Cint
    try
        for filename ∈ ARGS
            prove(filename)
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

end
