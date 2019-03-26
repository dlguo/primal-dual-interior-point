using MatrixDepot
using Test
using Printf
using SparseArrays
using LinearAlgebra
using Statistics

include("solver3.jl")
include("alpha_max.jl")
include("starting.jl")
include("convert2standard.jl")
include("presolve.jl")
include("solve_standardlp.jl")
include("phaseone.jl")

struct IplpSolution
    x::Vector{Float64} # the solution vector 
    flag::Bool         # a true/false flag indicating convergence or not
    cs::Vector{Float64} # the objective vector in standard form
    As::SparseMatrixCSC{Float64} # the constraint matrix in standard form
    bs::Vector{Float64} # the right hand side (b) in standard form
    xs::Vector{Float64} # the solution in standard form
    lam::Vector{Float64} # the solution lambda in standard form
    s::Vector{Float64} # the solution s in standard form
end  

struct IplpProblem
    c::Vector{Float64}
    A::SparseMatrixCSC{Float64} 
    b::Vector{Float64}
    lo::Vector{Float64}
    hi::Vector{Float64}
end

function convert_matrixdepot(mmmeta)
    return IplpProblem(
        vec(mmmeta.c),
        mmmeta.A,
        vec(mmmeta.b),
        vec(mmmeta.lo),
        vec(mmmeta.hi))
end

""" 
soln = iplp(Problem,tol) solves the linear program:

minimize c'*x where Ax = b and lo <= x <= hi

where the variables are stored in the following struct:

Problem.A
Problem.c
Problem.b   
Problem.lo
Problem.hi

and the IplpSolution contains fields 

[x,flag,cs,As,bs,xs,lam,s]

which are interpreted as   
a flag indicating whether or not the
solution succeeded (flag = true => success and flag = false => failure),

along with the solution for the problem converted to standard form (xs):

minimize cs'*xs where As*xs = bs and 0 <= xs

and the associated Lagrange multipliers (lam, s).

This solves the problem up to 
the duality measure (xs'*s)/n <= tol and the normalized residual
norm([As'*lam + s - cs; As*xs - bs; xs.*s])/norm([bs;cs]) <= tol
and fails if this takes more than maxit iterations.
"""

function iplp(Problem, tol; maxit=100)
    ### test input data
    
    @show m0,n0 = size(Problem.A)
    
    if length(Problem.b) != m0 || length(Problem.c) != n0 || length(Problem.lo) != n0 || length(Problem.hi) != n0
        DimensionMismatch("Dimension of matrices A, b, c mismatch. Check your input.")
    end

    @printf("Problem size: %d, %d\n",m0,n0)

    ### presolve stage

    Ps, ind0c, dup_main_c, ind_dup_c = presolve(Problem)

    ### convert to standard form
    @show size(Ps.A)
    @show rank(Array{Float64}(Ps.A))
    
    A,b,c,ind1,ind2,ind3,ind4 = convert2standard(Ps)
    @show size(A)
    @show rank(Array{Float64}(A))
    ### detect infeasibility

    if phaseone(A,b)
        @warn "This problem is infeasible."
        return IplpSolution(vec([0.]),false,vec(c),A,vec(b),vec([0.]),vec([0.]),vec([0.]))
    end

    @printf("\n=============== MPCIP solver ===============\n%3s %6s %11s %9s %9s\n", "ITER", "MU", "RESIDUAL", "ALPHAX", "ALPHAS")

    ### solve the original problem

    x1,lambda1,s1,flag,iter = solve_standardlp(A,b,c,maxit,tol,true)

    @printf("============================================\n")

    # @show iter

    x = get_x(Ps,ind1,ind2,ind3,ind4,x1)

    x = revProb(Problem, ind0c, dup_main_c, ind_dup_c, x)

    if flag == true
        @printf("This problem is solved with optimal value of %.2f.\n\n", dot(Problem.c, x))
    else
        @printf("\nThis problem does not converge in %d steps.", maxit)
    end

    return IplpSolution(vec(x),flag,vec(c),A,vec(b),vec(x1),vec(lambda1),vec(s1))
end
