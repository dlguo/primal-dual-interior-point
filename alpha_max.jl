"""
compute 
    arg max{ alpha in [0,hi] | x + alpha*dx >=0 }

"""
function alpha_max(x,dx,hi)
    n = length(x)

    alpha = -1.0

    for i=1:n
        if dx[i] < 0

            a = -x[i]/dx[i]

            if alpha < 0
                alpha = a
            else
                alpha = min(alpha,a)
            end
        end
    end

    if alpha < 0
        alpha = Inf
    end

    alpha = min(alpha,hi)

    return alpha
end
