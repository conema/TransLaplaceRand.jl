module TransLaplaceRand

using Roots
import Roots: bracket, adjust_bracket, NullTracks, isissue, log_message, AbstractUnivariateZeroState
import Roots.Accessors: @reset
using InverseLaplace: Fourier_series
using Random
using Distributions

include("newton_bracketed.jl")
include("rlaptrans.jl")

export rlaptrans

end
