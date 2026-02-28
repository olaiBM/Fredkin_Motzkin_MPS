using ITensors, LinearAlgebra
using ITensorMPS







function MotzkinMPS(N::Int, sites::Vector, q::Float64)
    """
    Returns the exact MPS representation of the Motzkin GS as an ITensor MPS

    N::Int                      # Length of Motzkin chain
    sites::Vector               # Vector containg the site indices
    q::Float64                  # Deformation parameter
    """

    @assert N%2 == 0 "Only works for even N"

    tensors = [ITensor() for _ in 1:N+2]
    internal_indices= create_internal_index_list_M(N)  

    tensors[1] = boundary_tensor_M(internal_indices[1], "start")
    tensors[end] = boundary_tensor_M(internal_indices[end], "end")

    for pos in 1:N
        tensors[pos+1] = A_M(sites[pos], internal_indices[pos], internal_indices[pos+1], N, q, pos)
    end

    tensors[2] = tensors[1] * tensors[2]
    tensors[end-1] = tensors[end] * tensors[end-1]
    tensors = tensors[2:end-1]
    return MPS(tensors)


end




function boundary_tensor_M(k::Index, position::String)
    """
    Returns boundary tensor, used to enforce fixed boundary conditions
    """
    delta = ITensor()

    if position == "start"
        delta = ITensor(k)
    elseif position == "end"
        delta = ITensor(dag(k))
    else
        error("Position must be even or odd")
    end

    delta[k => 1] = 1.0
    return delta
end

function create_internal_index_list_M(N::Int)
    
    """
    Creates list of internal_indices used to create the Motzkin MPS tensors (and boundary tensors)
    """
    internal_indices = [Index(QN("init", 1) => 1;tags = "init link $i") for i in 1:N+1]
    internal_indices[1] = Index(QN("Sz", 0) => 1; tags= "h link = 0")
    internal_indices[N+1] = Index(QN("Sz", 0) => 1; tags= "h link = $N")
    for pos in 1:(N-1)
        max_QN = min(pos, N - pos)
        start_QN = 0
        internal_indices[pos + 1] = Index([QN(("Sz", qn)) => 1 for qn in start_QN:max_QN]; tags="h link=$pos")
    end
    return internal_indices

end




function A_M(site_ind::Index, i_m::Index, i::Index, N::Int, q::Float64, pos::Int)
    """
    Creates the Motzkin MPS tensor
    """
    A = ITensor() # Tensor to be returned

    A = ITensor(dag(i_m), dag(site_ind), i)

    min_QN = 0   # is zero 
    max_QN = min(pos-1, N - pos+1) # Max QN for the i_m index
    for qn in min_QN:(Int(max_QN))

        if pos <= N/2
            if qn == 0
                ### Then we have either a flat spin or up spin (can go into negative qn sector)

                A[i_m => (qn+1), site_ind => 1, i => (qn+2)] = tensor_val_M("up", q, N, qn, pos) 
                A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
            else
                ### Then we can go either of the three steps
                A[i_m => (qn+1), site_ind => 1, i => (qn+2)] = tensor_val_M("up", q, N, qn, pos) 
                A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
                A[i_m => (qn+1), site_ind => 3, i => (qn)] = tensor_val_M("down", q, N, qn, pos) 

            end
        else

            if qn == max_QN
                ### Then we can only go down
                A[i_m => (qn+1), site_ind => 3, i => (qn)] = tensor_val_M("down", q, N, qn, pos) 
            elseif qn == (max_QN-1)
                ### Then we can only go flat or down
                if qn == 0
                    A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
                else
                    A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
                    A[i_m => (qn+1), site_ind => 3, i => qn] = tensor_val_M("down", q, N, qn, pos) 
                end
            else
                if qn == 0
                    ### Then we have either a flat spin or up spin (can go into negative qn sector)
                    
                    A[i_m => (qn+1), site_ind => 1, i => (qn+2)] = tensor_val_M("up", q, N, qn, pos) 
                    A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
                else
                
                    A[i_m => (qn+1), site_ind => 1, i => (qn+2)] = tensor_val_M("up", q, N, qn, pos) 
                    A[i_m => (qn+1), site_ind => 2, i => (qn+1)] = tensor_val_M("flat", q, N, qn, pos) 
                    A[i_m => (qn+1), site_ind => 3, i => (qn)] = tensor_val_M("down", q, N, qn, pos) 
                end

            end
        end

    end    
    return A
end



function tensor_val_M(spin::String, q::Float64, N::Int, qn::Int, pos::Int)
    """
    Provides the value of the tensor entries, using that we consider the minimum (resp. maximum) area 
    walk as the walk of unit area for q<1 (resp. q>1).
    """
    entry = 0                    # entry of tensor to be returned
    exponent = 0
    max_QN_im = min(pos-1, N - pos +1)
    max_QN_i = min(pos, N - pos)
    if q> 1.0
        ### We then consider the maximum area the walk of unit area (so other walks are suppressed)
        if spin == "up"
            if pos <= N/2
                exponent = 2*(qn - max_QN_im) 
            else
                exponent = 2*(qn+1 - max_QN_i) - 2
            end
        elseif spin == "flat"
            if pos <= N/2
                
                exponent = 2*(qn - max_QN_im) - 1
            else
                exponent = 2*(qn- max_QN_i) - 1
            end

        elseif spin == "down"
            if pos <= N/2
                exponent = 2*(qn - max_QN_im) - 2
            else
                exponent = 2*(qn - max_QN_im)
            end
        else
            error("Spin can only be up, flat or down")
        end
        entry = q^(exponent)
    elseif q == 1
        ### Then all walks are weighted equally
        entry = 1.0
    else
        ### We then consider the minimum area the walk of unit area (so other walks are suppressed)
        if spin == "up"
            
            exponent = 2*qn + 1
        elseif spin == "flat"
            
            exponent = 2*qn 
        elseif spin == "down"
            exponent =2*(qn-1) +1
        else
            error("Spin can only be up, flat or down")
        end
        entry = q^(exponent)


    end

    return entry
end


