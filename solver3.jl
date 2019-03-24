# unreduced form

function fact3(A,x,s)
    m,n = size(A)

    M = [spzeros(m,m) A spzeros(m,n);
     A' spzeros(n,n) speye(n);
      spzeros(n,m) spdiagm(s[:,1]) spdiagm(x[:,1])]

    f = lufact(M)

    return f
end

function solve3(f,A,x,s,rb,rc,rxs)
    m = length(rb)
    n = length(rc)

    b = full([-rb; -rc; -rxs])
    b = f\b

    dlam = b[1:m]
    dx = b[1+m:m+n]
    ds = b[1+m+n:m+2*n]

    return dlam,dx,ds
end