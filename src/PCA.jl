using LinearAlgebra

function PCA(filename::AbstractString)

	q = 21
	Pi, Pij, Pi_tens, Pij_tens, C = compute_freq_corr(filename)

	v = eigvecs(C)
	l = eigvals(C)

	# max eigenvec/val 
	v1 = v[:,end]
	l1 = l[end]
	# second max
	v2 = v[:,end-1]
	l2 = l[end-1]
	# third max
	v3 = v[:,end-2]
	l3 = l[end-2]

	#make Z as 0 and 1
	max_gap_fraction = 0.9
	Z = GaussDCA.ReadFastaAlignment.read_fasta_alignment(filename, max_gap_fraction)
	N,M = size(Z)
	Z_tmp = zeros(N*q, M)
	for m = 1:M
		for n = 1:N
			Z_tmp[Z[n,m] + q*(n-1), m] = 1
		end
	end
	
	#compute PCA componenets
	x1 = Z_tmp' * v1
	x2 = Z_tmp' * v2
	x3 = Z_tmp' * v3
	
	return x1,x2,x3
end
