"""
conver the lp problem
min c^T*x
s.t. Ax = b
lo <= x <= hi

to an lp in standard form
min c^T*x
s.t. Ax = b
x >= 0
"""
function convert2standard(P)
    inf = 1.0e300

    m,n0 = size(P.A)

    index1 = zeros(Int64,0)
    index2 = zeros(Int64,0)
    index3 = zeros(Int64,0)
    index4 = zeros(Int64,0)
    n = zeros(Int64,4)

    #=    lo   hi
    n1 :       
    n2 :   x 
    n3 :        x
    n3 :   x    x
    =#

    for i=1:n0
        if P.lo[i] < -inf
            if P.hi[i] > inf
                n[1] += 1
                index1 = [index1;i]
            else
                n[3] += 1
                index3 = [index3;i]
            end
        else
            if P.hi[i] > inf
                n[2] += 1
                index2 = [index2;i]
            else
                n[4] += 1
                index4 = [index4;i]
            end
        end
    end

    cs = [P.c[index1,1]; -P.c[index1,1]; P.c[index2,1];
     -P.c[index3,1]; P.c[index4,1]; zeros(n[4],1)]

    As = [P.A[:,index1] -P.A[:,index1] P.A[:,index2] -P.A[:,index3] P.A[:,index4] spzeros(m,n[4]);
     spzeros(n[4],2*n[1]+n[2]+n[3]) Matrix{Float64}(I,n[4],n[4]) Matrix{Float64}(I,n[4],n[4])]

    bs = [P.b-P.A[:,index2]*P.lo[index2,1]-P.A[:,index3]*P.hi[index3,1]-P.A[:,index4]*P.lo[index4,1];
     P.hi[index4,1]-P.lo[index4,1]]

    return As,bs,cs,index1,index2,index3,index4
end

"""
obtain the solution vector of the original problem
"""

function get_x(P,ind1,ind2,ind3,ind4,xs)
    m,n = size(P.A)
    n1 = length(ind1)
    n2 = length(ind2)
    n3 = length(ind3)
    n4 = length(ind4)

    x = zeros(n)

    x[ind1] = xs[1:n1]-xs[1+n1:2*n1]
    x[ind2] = xs[1+2*n1:n2+2*n1]+P.lo[ind2]
    x[ind3] = P.hi[ind3]-xs[1+2*n1+n2:n3+2*n1+n2]
    x[ind4] = xs[1+2*n1+n2+n3:n4+2*n1+n2+n3]+P.lo[ind4]

    return x
end
