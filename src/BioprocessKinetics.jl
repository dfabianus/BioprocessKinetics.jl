module BioprocessKinetics

# Write your package code here.
export monod, blackman, teissier, moser

"Only the substrate uptake step is rate limiting"
function monod(c::AbstractFloat, k::AbstractFloat)
    c / (k + c)
end

"Another rate limiting step besides substrate uptake determines the maximum rate"
function blackman(c::AbstractFloat, k::AbstractFloat)
    minimum([1, k * c])
end

"Empirical equation"
function teissier(c::AbstractFloat, k::AbstractFloat)
    1 - exp(-k * c)
end

"Substrate uptake with higher order of reaction, e.g. for gaseous substrates"
function moser(c::AbstractFloat, k::AbstractFloat, N::Int)
    c^N / (k^N + c^N)
end

"Diffusion layer around the cell"
function contois(cₛ::AbstractFloat, cₓ::AbstractFloat, kₛ::AbstractFloat)
    cₛ / (kₛ * cₓ + cₛ)
end

"Considers back diffusion of inner substrate"
function powel(c::AbstractFloat, kₛ::AbstractFloat, k₁::AbstractFloat, r₁::AbstractFloat)
    (c - k₁ * r₁) / (k + c - k₁ * r₁)
end

function haldane(c::AbstractFloat, kₛ::AbstractFloat, kᵢ::AbstractFloat)
    c / (kₛ + c * (1 + c / kᵢ))
end
end
