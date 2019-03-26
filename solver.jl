# augmented system

"""
sparse LDL factorization of 
| -X^(-1)*S  A^T |
|     A       0  |

returns the factors given by ldltfact
"""
function fact(A,x,s)
    m,n  = size(A)

    l1 = 1.0e16
    l2 = 1.0e-15

    d = s./x
    d = min(d,l1)
    # d = max(d,l2)
    d = -d

    SID = [ spzeros(m,m) A; A' sparse(collect(1:n),collect(1:n),d[:,1])]

    f = ldlt(SID)

    # try f = ldltfact(SID)
    # catch
    #     MatrixMarket.mmwrite("sid.mtx",SID)
    #     exit
    # end

    return f
end


"""
This function solves the following linear system

|  0  A  0 |   | dlam |   | -r_b  |
| A^T 0  I | * |  dx  | = | -r_c  |
|  0  S  X |   |  ds  |   | -r_xs |

by solving 
| -X^(-1)*S  A^T | * |  dx  | = | -r_c+X^(-1)*r_xs |
|     A       0  |   | dlam |   |       -r_b       |

and 

ds = -X^(-1)*r_xs - X^(-1)*S*dx

it returns dlam,dx,ds
"""

function solve(f,x,s,rb,rc,rxs)
    n = length(x)
    m = length(rb)

    # b1 = [-rb; -rc+rxs./x]
    b1 = Array{Float64}([-rb; -rc+rxs./x])

    x1 = f\b1

    dlam = x1[1:m]
    dx   = x1[1+m:m+n]

    ds = -(rxs+s.*dx)./x

    return dlam,dx,ds
end
