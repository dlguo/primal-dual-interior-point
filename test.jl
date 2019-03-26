#!/usr/bin/julia

# using MatrixDepot
# using Base.Test

include("iplp.jl")
# include("solver.jl")
# include("solver2.jl")
# include("solver3.jl")
# include("alpha_max.jl")
# include("starting.jl")
# include("convert2standard.jl")
# include("presolve.jl")
# include("solve_standardlp.jl")
# include("phaseone.jl")

problem_list = ["lp_afiro","lp_brandy","lp_fit1d","lp_adlittle",
"lp_agg","lp_ganges","lp_stocfor1", "lp_25fv47", "lpi_chemcom"]

for i = 1:length(problem_list)
    @printf("%d. %s\n", i, problem_list[i])
end
@printf("0. other\n")
@printf("Which problem to solve? ")

k = parse(Int,readline(stdin))

if k == 0
    @printf("Please enter the problem name (e.g. lp_afiro): ")
    name = "LPnetlib/"*readline(stdin)[1:end-1]
else
    name = "LPnetlib/"*problem_list[k]
end

@printf("Solving %s.", name)

try global P = convert_matrixdepot(mdopen(name))
catch 
    mdopen(name)
    global P = convert_matrixdepot(mdopen(name))
end

tol=1e-8
solution = @time iplp(P, tol; maxit=100)


# @show solution.flag
# @show dot(P.c,solution.x)
