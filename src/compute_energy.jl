#The dimension of the input parameters must be h=(q,N) and J=(q,q,N,N)

#compute energy of a single sequence
#h=(q,N) and J=(q,N)
function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)

    N = size(h)[2]
    q = size(h)[1]
    E = 0.0
    for i = 1:N
        E -= h[S[i],i]
        for j = (i+1):N 
			E -= J[S[i],S[j],i,j]
		end
	end
return E
end

#compute energy of a MSA
#S is a matrix M x N (M=num of sequence, N=num of spins)
function compute_energy_MSA(h::Array{Float64,2},
                            J::Array{Float64,4},
                            S::Matrix)

    M,N = size(S)
    final_E = SharedArray{Float64}(M) 
    q=size(h)[1]

    println("N=", N, " M=",M, " q=",q)

    @sync begin
        @parallel for i_seq = 1:M
            println("M ",i_seq)
            final_E[i_seq] = compute_energy_single_sequence(h,J,S[i_seq,:])
        end
    end


    return final_E
end

#compute energy of a MSA
#fastafile is a fasta file 
function compute_energy_MSA(h::Array{Float64,2},
                            J::Array{Float64,4},
                            fastafile::AbstractString)

    S = fasta2matrix(fastafile,0.9)
    return compute_energy_MSA(h,J,S)
end

