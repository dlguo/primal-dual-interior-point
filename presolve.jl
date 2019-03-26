function presolve(P)
"""
Presolve function

"""
    m, n = size(P.A)
    ind0r = zeros(Int64,0)
    ind0c = zeros(Int64,0)
    ind_dup_r = zeros(Int64,0)
    ind_dup_c = zeros(Int64,0)
    dup_main_c = Array[]
    dup_del_c = Array[]
    remove_ind_row = zeros(Int64,0)
    remove_ind_col = zeros(Int64,0)

    # zero rows and columns in A
    # zero rows
    for i = 1:m
        j = 1
        while (j <= n) && (P.A[i,j] == 0.0)
            j += 1
        end

        if j == n+1
            if P.b[i] == 0.0
                ind0r = [ind0r; i]
            else
                warn("This problem is infeasible.")
                return false
            end
        end
    end
    # zero columns
    for j =1:n
        i = 1
        while (i <= m) && (P.A[i,j] == 0.0)
            i += 1
        end

        if i == m+1
            ind0c = [ind0c; j]
        end
    end

    # duplicated rows
    for i=1:(m-1)
        if (i in ind_dup_r) || (i in ind0r)
            continue
        end

        for j=(i+1):m
            if (j in ind_dup_r) || (j in ind0r)
                continue
            end

            k = 1
            while (k <= n) && (P.A[i,k] == P.A[j,k])
                k += 1
            end

            if k == n+1
                if P.b[i] == P.b[j]
                    ind_dup_r = [ind_dup_r; j]
                else
                    warn("This problem is infeasible.")
                    return false
                end
            end

        end
    end
    @show remove_ind_row = [ind0r; ind_dup_r]

    # duplicate columns 
    ### need modification
    for i=1:(n-1)
        dup_item = Array(i:i)
        for j=(i+1):n
            if !(j in ind_dup_c) && !(j in ind0c) && (P.A[:, i] == P.A[:, j])
                dup_item = [dup_item; j]
            end
        end
        if length(dup_item)>1
            minv, mini = findmin(P.c[dup_item])
            dup_flags = trues(length(dup_item))
            dup_flags[mini] = false
            dup_main_c = [dup_main_c; dup_item[mini]]
            ind_dup_c = [ind_dup_c; dup_item[dup_flags]]
        end
    end
    @show remove_ind_col = [ind0c; ind_dup_c]

    indpre_row = Array(1:m)
    indpre_col = Array(1:n)
    flags_row = trues(m)
    flags_col = trues(n)
    flags_row[remove_ind_row] .= false
    flags_col[remove_ind_col] .= false

    return IplpProblem(P.c[flags_col], P.A[flags_row, flags_col],
     P.b[flags_row], P.lo[flags_col], P.hi[flags_col]), ind0c, dup_main_c, ind_dup_c
end

function revProb(P, ind0c, dup_main_c, ind_dup_c, x1)
    m, n = size(P.A)
    x = Array{Float64}(undef, n)
    fill!(x, Inf)
    j = 1
    for i = 1:n
        if x[i] == Inf
            if i in ind0c
                if P.c[i] > 0
                    x[i] = P.lo[i]
                elseif P.c[i] <0
                    x[i] = P.hi[i]
                else
                    x[i] = 0.
                end
            elseif i in dup_main_c
                x[i] = x1[j]
                j += 1
            elseif i in ind_dup_c
                x[i] = 0.
            else
                x[i] = x1[j]
                j += 1
            end
        end
    end
    @test j == length(x1) + 1
    
    return x
end
