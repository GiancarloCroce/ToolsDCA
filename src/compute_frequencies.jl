using FastaIO
using GaussDCA

function compute_freq_corr(fastafile::AbstractString; pseudocount::Real=0.8)
	max_gap_fraction = 0.9
	theta = 0.0 #standard is 0.2
	q = 21
	#Z = ToolsDCA.fasta2matrix(fastafile,max_gap_fraction) 
	Z = GaussDCA.read_fasta_alignment(fastafile, max_gap_fraction)
	N,M = size(Z)
	Pi_true, Pij_true, Meff, _ = compute_new_frequencies(Z, q+1, theta)
	#q+1 to get Pi = q*N (not (q-1)*N)
	Pi, Pij = GaussDCA.add_pseudocount(Pi_true, Pij_true, Float64(pseudocount), N, q+1)
	#make tensor form 
	Pi_tens = reshape(Pi,q,N)
	Pij_tens = zeros(q,q,N,N)
	for i  = 1:N
		for j = (i+1):N
			Pij_tens[:,:,i,j] = Pij[q*(i-1)+1 : q*(i), q*(j-1)+1 : q*(j)]
			Pij_tens[:,:,j,i] = Pij[q*(j-1)+1 : q*(j),q*(i-1)+1 : q*(i)]
		end
	end
	C = GaussDCA.compute_C(Pi, Pij)
	return Pi,Pij,Pi_tens,Pij_tens,C
end


function compute_new_frequencies(Z::Matrix{Int8}, q, theta)
    W, Meff = GaussDCA.compute_weights(Z, q, theta)
    Pi_true, Pij_true = compute_freqs(Z, W, Meff,q)
    return Pi_true, Pij_true, Meff, W
end


function compute_freqs(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64, q::Int64)
    N, M = size(Z)
    s = q - 1
    Ns = N * s
    Pij = zeros(Ns, Ns)
    Pi = zeros(Ns)
    ZZ = Vector{Int8}[vec(Z[i,:]) for i = 1:N]
    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        for k = 1:M
            a = Zi[k]
            a == q && continue
            Pi[i0 + a] += W[k]
        end
        i0 += s
    end
    Pi /= Meff
    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        j0 = i0
        for j = i:N
            Zj = ZZ[j]
            for k = 1:M
                a = Zi[k]
                b = Zj[k]
                (a == q || b == q) && continue
                Pij[i0+a, j0+b] += W[k]
            end
            j0 += s
        end
        i0 += s
    end
    for i = 1:Ns
        Pij[i,i] /= Meff
        for j = i+1:Ns
            Pij[i,j] /= Meff
            Pij[j,i] = Pij[i,j]
        end
    end
    return Pi, Pij
end

function add_pseudocount(Pi_true::Vector{Float64}, Pij_true::Matrix{Float64}, pc::Float64, N::Int, q::Int)
    pcq = pc / q
    Pij = (1 - pc) * Pij_true .+ pcq / q
    Pi = (1 - pc) * Pi_true .+ pcq
    s = q - 1
    i0 = 0
    for i = 1:N
        xr::UnitRange{Int} = VERSION < v"0.7.0-DEV.1759" ? i0 + (1:s) : i0 .+ (1:s)
	Pij[xr, xr] = (1 - pc) * Pij_true[xr, xr]
        for alpha = 1:s
            x = i0 + alpha
            Pij[x, x] += pcq
	end
        i0 += s
    end
    return Pi, Pij
end

compute_C(Pi::Vector{Float64}, Pij::Matrix{Float64}) = Pij - Pi * Pi'


