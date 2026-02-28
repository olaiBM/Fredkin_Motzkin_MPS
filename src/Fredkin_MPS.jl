using ITensors, LinearAlgebra
using ITensorMPS



function FredkinMPS(N::Int, sites::Vector, q::Float64)
    """
    Returns the exact MPS representation of the Fredkin GS as an ITensor MPS

    N::Int                      # Length of Fredkin chain
    sites::Vector               # Vector containg the site indices
    q::Float64                  # Deformation parameter
    """
    
    @assert N % 2 == 0 "N must be even for the Fredkin chain"
   
    tensors = [ITensor() for _ in 1:N+2]
    internal_indices = create_internal_index_list_F(N)  

    tensors[1] = boundary_tensor_F(internal_indices[1],"start")
    tensors[end] = boundary_tensor_F(internal_indices[end], "end")

    for pos in 1:N
        if pos % 2 == 1
            tensors[pos+1] = A_o_F(sites[pos], internal_indices[pos], internal_indices[pos+1], N, q, pos)
        else
            tensors[pos+1] = A_e_F(sites[pos], internal_indices[pos], internal_indices[pos+1], N, q, pos)
        end
    end
    tensors[2] = tensors[1] * tensors[2]
    tensors[end-1] = tensors[end] * tensors[end-1]
    tensors = tensors[2:end-1] # Return the MPS tensors (not boundary tensors at 1 and)
    return MPS(tensors)

    
end



function boundary_tensor_F(k::Index, position::String)
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

function create_internal_index_list_F(N::Int)
    """
    Creates list of the internal indices used to create the Fredkin MPS tensors (and boundary tensors)
    """
    internal_indices = [Index(QN("init", 1) => 1;tags = "init link $i") for i in 1:N+1]
    internal_indices[1] = Index(QN("Sz", 0) => 1; tags= "h link = 0")
    internal_indices[N+1] = Index(QN("Sz", 0) => 1; tags= "h link = $N")
    for pos in 1:(N-1)
        max_QN = min(pos, N - pos)
        start_QN = pos%2
        internal_indices[pos + 1] = Index([QN(("Sz", qn)) => 1 for qn in start_QN:2:max_QN]; tags="h link=$pos")
    end
    return internal_indices

end




function A_o_F(site_ind::Index, i_m::Index, i::Index, N::Int, q::Float64, pos::Int)
    """
    Creates the Fredkin MPS tensor appropriate for odd sites
    """
    A = ITensor() # Tensor to be returned

    A = ITensor(dag(i_m), dag(site_ind), i)

    start_QN = 0   # For odd sites, the i_m index starts at 0
    max_QN = min(pos-1, N - pos+1) # Max QN for the i_m index
    index = 1  

    for qns in start_QN:2:(Int(max_QN))

        if pos <= N/2
            if index == 1
                A[i_m => index, site_ind => 1, i => index] = tensor_val_F("up", q, N, qns, pos) 
            else
                A[i_m => index, site_ind => 1, i => index] = tensor_val_F("up", q, N, qns, pos) 
                A[i_m => index, site_ind => 2, i => index-1] = tensor_val_F("down", q, N, qns-1, pos) 
            end

        else 
            if index == 1
                A[i_m => index, site_ind => 1, i => index] = tensor_val_F("up", q, N, qns, pos) 
            elseif (index -1) * 2 == max_QN   # Have reached max QN for i_m, must take down step when in second half of spin chain
                A[i_m => index, site_ind => 2, i => index-1] = tensor_val_F("down", q, N, qns-1, pos)
            else
                A[i_m => index, site_ind => 1, i => index] = tensor_val_F("up", q, N, qns, pos) 
                A[i_m => index, site_ind => 2, i => index-1] = tensor_val_F("down", q, N, qns-1, pos) 

            end

        end
        
        index += 1
    end    
    return A
end


function A_e_F(site_ind::Index, i_m::Index, i::Index, N::Int, q::Float64, pos::Int)
    """
    Creates the Fredkin MPS tensor appropriate for even sites
    """
    A = ITensor()

    A = ITensor(dag(i_m), dag(site_ind), i)

    start_QN = 1 # Is always one
    max_QN = min(pos-1, N - pos+1)
    index = 1

    for qns in start_QN:2:(Int(max_QN))

        if pos <= N/2

            A[i_m => index, site_ind => 1, i => index + 1] = tensor_val_F("up", q, N, qns, pos)# Spin up part
            A[i_m => index, site_ind => 2, i => index] = tensor_val_F("down", q, N,  qns-1, pos) # Spin down part

        else
            if ((index -1)*2 + 1) == max_QN
                A[i_m => index, site_ind => 2, i => index] = tensor_val_F("down", q, N, qns-1, pos) # Spin down part
            else
                A[i_m => index, site_ind => 1, i => index + 1] = tensor_val_F("up", q, N, qns, pos) # Spin up part
                A[i_m => index, site_ind => 2, i => index] = tensor_val_F("down", q, N, qns-1, pos) # Spin down part

            end

        end
        
        index += 1
    end    
    return A
end


function exponent_qg1(spin::String, pos::Int, N::Int)
    """
    Returns the exponent to assure "normalized state" for the q>1 case
    """

    exponent = 0
    if pos <= Int(N/2)
        if spin == "down"
            N_squares = (N/2) - pos +0.5
            exponent = -2*N_squares 
        end
    else
        if spin == "up"
            N_squares = pos - (N/2) - 1 +0.5
            exponent = -2*N_squares
        end
    end
    return exponent
end

function tensor_val_F(spin::String, q::Float64, N::Int, qn::Int, pos::Int)
    """
    Returns the value of the Fredkin MPS tensor entry for a given index configuration (qn)
    """
    component= 0 # component to be returned
    if q> 1.0
        exponent = exponent_qg1(spin, pos, N) 
        component = q^(exponent)
    elseif q == 1
        component = 1.0
    else

        exponent = qn
        component = q^(exponent)


    end

    return component
end



