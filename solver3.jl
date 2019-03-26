# unreduced form

function fact3(A,x,s)
    m,n = size(A)

    M = [spzeros(m,m) A spzeros(m,n);
     A' spzeros(n,n) Matrix{Float64}(I,n,n);
      spzeros(n,m) spdiagm(0=>s[:,1]) spdiagm(0=>x[:,1])]

    f = lu(M)

    return f
end

function solve3(f,A,x,s,rb,rc,rxs)
    m = length(rb)
    n = length(rc)

    b = Array{Float64}([-rb; -rc; -rxs])
    b = f\b

    dlam = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:m+2*n]

    return dlam,dx,ds
end