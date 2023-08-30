using TransLaplaceRand: rlaptrans
using Distributions
using FLoops

function generate_z(draws, lambda)
    z_pdf(t) = exp(-t^lambda)

    rlaptrans(draws, z_pdf)
end

function calculate_prob_product_max(draws, zs, utilities, lambdas, nest, product)
    levels = [0]

    for n in utilities
        push!(levels, length(n))
    end

    @floop for draw in 1:draws
        Us = []

        for n in eachindex(lambdas)
            gumbels = rand(Gumbel(), length(utilities[n]))

            e = lambdas[n] .* (log(zs[n, draw]) .+ gumbels )
            U = utilities[n] .+ e

            append!(Us, U)
        end

        max_u = maximum(vcat(Us...))

        if max_u == Us[levels[nest] + product]
            @reduce is_max += 1
        end
    end

    return is_max / draws
end

function generate_z_nests(draws, lambdas)
    generated_lambdas = zeros(length(lambdas), draws)

    @floop for nest_i in eachindex(lambdas)
        @inbounds generated_lambdas[nest_i, :] = generate_z(draws, lambdas[nest_i])
    end

    return generated_lambdas
end

function calculate_prob_product_closed(utilities, lambdas, nest, product)
    in_nest(nest, utilities, lambdas) = sum(exp.(utilities[nest] ./ lambdas[nest]))

    all_nest(nests, utilities, lambdas) = sum(
        [in_nest(nest, utilities, lambdas) ^ lambdas[nest] for nest in nests]
    )

    in_nest_power = in_nest(nest, utilities, lambdas) ^ (lambdas[nest] - 1)

    all_nest_res = all_nest(1:length(utilities), utilities, lambdas)

    return (exp(utilities[nest][product] / lambdas[nest]) * in_nest_power) / all_nest_res
end

draws = 10000
lambdas = [0.8, 0.75, 0.2]
utilities = [
    [5, 1, 3],
    [2, 4, 6],
    [1, 2]
]

nest = 1
product = 1

zs = generate_z_nests(draws, lambdas)

println(
    calculate_prob_product_closed(utilities, lambdas, nest, product)
)

println(
    calculate_prob_product_max(draws, zs, utilities, lambdas, nest, product) 
)