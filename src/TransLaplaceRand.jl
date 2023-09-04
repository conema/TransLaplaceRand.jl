module TransLaplaceRand

using Roots
import Roots: @set!, bracket, adjust_bracket, NullTracks, isissue, log_message, AbstractUnivariateZeroState
using InverseLaplace: Fourier_series
using Random
using Distributions

include("newton_bracketed.jl")
include("rlaptrans.jl")

export rlaptrans

end
