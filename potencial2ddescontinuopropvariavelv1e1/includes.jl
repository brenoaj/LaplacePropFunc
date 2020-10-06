#=
cd C:/Users/Admin/Documents/julia/julia/
=#
using Plots
using LinearAlgebra
using FastGaussQuadrature
using ForwardDiff
using TimerOutputs
using Statistics
using FFTW
# using StatsPlots
#using Triangle
# using WriteVTK


include("dad.jl") # Arquivos de dados de problemas gerais (dad_0, dad_1, ..., dad_9)
include("format.jl")
include("gera_p_in.jl")
include("inpoly.jl")
include("telles.jl")
include("cal_HeG.jl")
include("calc_fforma.jl")
include("mecmatrizvetor.jl")
include("pontoiinterno.jl")
include("Contorno.jl")
include("transforma.jl")
include("mostra_problema.jl")

mutable struct dados
    NOS ::Array{Float64,2}
    ELEM ::Array{UInt32,2}
    tipoCDC ::Array{Float64,1} #Why is this float?
    valorCDC ::Array{Float64,1}
    afasta ::Float64
    k ::Float64
    kf ::Function
end
