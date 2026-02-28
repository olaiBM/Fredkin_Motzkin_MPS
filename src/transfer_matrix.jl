using LinearAlgebra



function tridiag_transfer_matrix_F(L::Int, q::Float64)
    """
    Returns the maximum dimensional block of the tridiagonal Fredkin transfer matrix
    """
    d = zeros(L+1)  # Diagonal 
    dl = zeros(L)  # Subdiagonal

    exponent = 0
    weight = 0
    for i in 1:Int(L)
        if q>1
            exponent = (2 * i) - 1 - (2 * L) 
        else
            exponent = (2 * i) - 1
        end
        weight = q^exponent
        dl[i] = weight  # Only first N/2 values are nonzero
    end


    T = SymTridiagonal(d, dl)
    return T
end


function tridiag_transfer_matrix_M(L::Int, q::Float64)
    """
    Returns the maximum dimensional block of the tridiagonal Fredkin transfer matrix
    """
    d = zeros(L+1)  # Diagonal 
    dl = zeros(L)  # Subdiagonal

    exponent = 0
    weight = 0
    for i in 1:Int(L)
        if q>1
            exponent = (2 * i) - 1 - (2 * L) 
        else
            exponent = (2 * i) - 1
        end
        weight = q^exponent
        dl[i] = weight  # Only first N/2 values are nonzero
    end

    for i in 1:Int(L+1)
        if q>1
            exponent = 2 * (i - 1) - (2 * L)
        else
            exponent = 2 * (i - 1) 
        end
        weight = q^exponent
        d[i] = weight
    end
    # Construct symmetric tridiagonal matrix
    T = SymTridiagonal(d, dl)
    return T
end



function compute_gap(T::SymTridiagonal)
    """
    Computes the gap between the two largest eigenvalues
    """

    lambdas = sort(eigvals(T); rev=true)  # Descending order
    gap = (lambdas[1] - lambdas[2])/lambdas[1]
    return gap
end

function correlation_length(T::SymTridiagonal{Float64, Vector{Float64}})
    """
    Computes the correlation length from the transfer matrix
    """
    evals = eigvals(T)
    evals_sorted = sort(evals, rev=true)

    lambda_1, lambda_2 = evals_sorted[1:2] # The first one is doubly degenerae
    #println(lambda_1,"\t", lambda_2)
    return 1 / log(lambda_1 / abs(lambda_2))
end