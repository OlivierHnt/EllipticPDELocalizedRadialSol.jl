module EllipticPDELocalizedRadialSol

    import JLD2: load
    using RadiiPolynomial
    import RadiiPolynomial.LinearAlgebra: Diagonal, diag

include("klein_gordon.jl")
    export KleinGordon

include("swift_hohenberg.jl")
    export SwiftHohenberg

include("fitzhugh_nagumo.jl")
    export FitzHughNagumo

include("proof.jl")
    export prove

end
