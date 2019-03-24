# normal equation
# Cholesky factorization (skip small pivots)

function fact2(A,x,s)
    m,n = size(A)
    tol = 1.0e-12

    d = x./s
    L = A'

    for i = 1:n
        L[i,:] = d[i]*L[i,:]
    end

    L = A*L

    ind_skip = zeros(Int64,0)
    d_max = 0.0

    for i = 1:m

        if L[i,i] < tol*d_max
            ind_skip = [ind_skip; i]

            continue
        end

        d_max = max(d_max,L[i,i])

        L[i,i] = sqrt(L[i,i])

        if i == m
            break
        end

        for j = i+1:m
            L[i,j] = L[i,j]/L[i,i]

            for k = i+1:j
                L[k,j] = L[k,j] - L[i,j]*L[i,k]
            end
        end

    end

    # @show ind_skip

    return L,ind_skip
end

function solve2(L,ind_skip,A,x,s,rb,rc,rxs)
    m,n = size(A)

    dlam = (rxs-x.*rc)./s
    dlam = A*dlam-rb

    # back solve
    for i = 1:m
        if i in ind_skip
            dlam[i] = 0.0

            continue
        end

        dlam[i] = dlam[i]/L[i,i]

        if i == m 
            break
        end

        for j = i+1:m
            # dlam[i+1:m] = dlam[i+1:m] - L[i,i+1:m]'*dlam[i]
            dlam[j] = dlam[j] - L[i,j]*dlam[i]
        end
    end

    # forward solve
    for i = m:-1:1
        if i in ind_skip
            dlam[i] = 0.0

            continue
        end

        dlam[i] = dlam[i]/L[i,i]

        if i == 1
            break
        end

        for j = 1:i-1
            # dlam[1:i-1] = dlam[1:i-1] - L[1:i-1,i]*dlam[i]
            dlam[j] = dlam[j] - L[j,i]*dlam[i]
        end
    end

    ds = -rc - A'*dlam
    dx = -(rxs+x.*ds)./s

    return dlam,dx,ds
end

function invp(p)
    # inverse of a permutation
    n = length(p)

    q = zeros(Int64,n)

    for i = 1:n
        q[p[i]] = i
    end

    return q
end