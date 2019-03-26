"""
obtain the starting point for PD method
(page 410 on Wright)
"""

function starting_point(A,b,c)
    AA = A*A'

    f = cholesky(AA)
    # f = ldltfact(AA)
    # f = factorize(AA)

    # tilde
    x = f\b
    x = A'*x

    lambda = A*c
    lambda = f\lambda

    s = A'*lambda
    s = c-s

    # hat
    dx = max(-1.5*minimum(x),0.0)
    ds = max(-1.5*minimum(s),0.0)

    x = x.+dx
    s = s.+ds

    # ^0
    xs = dot(x,s)/2.0

    dx = xs/sum(s)
    ds = xs/sum(x)

    x = x.+dx
    s = s.+ds

    return x,lambda,s
end
