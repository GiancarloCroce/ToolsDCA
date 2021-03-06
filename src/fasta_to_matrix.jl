#it converts a fasta file to a matrix Z = MxN by mapping the letters (residues) to numbers

using Random 
using FastaIO
using Compat


#shuffle matrix 
function shuffle_matrix(a::Matrix,
                        dim::Int,
                        num_times::Int = 1)
    M,N = size(a)
    println("shuffling matrix with: N=",N," M=",M)
    b_old = a
    b_new = zeros(size(a))
    if dim == 1
        println("shuffle along dimension 0 (shuffling all elements of the matrix)")
        b = Random.shuffle(vec(a))
        return reshape(b,size(a))
    elseif dim == 1
        println("shuffle along dimension 1 (keep Pi, destroy Pij)")
        for ii = 1:num_times
            for k = 1:N
                b_new[:,k] = Random.shuffle(b_old[:,k])
            end
            b_old = b_new
        end
        return b_new
    elseif dim==2
        println("shuffle along dimension 2") 
        for ii = 1:num_times
            for k = 1:M
                b_new[k,:] = Random.shuffle(b_old[k,:])
            end
            b_old = b_new
        end
        return b_new
   end
end

#read a fasta and get a matrix MxN
function fasta2matrix(filename::AbstractString, max_gap_fraction::Real)
    f = FastaReader(filename)
    max_gap_fraction = Float64(max_gap_fraction)
    # pass 1
    seqs = Int[]
    inds = Int[]
    fseqlen = 0
    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end
    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")
    # pass 2
    Z = Array{Int8}(undef,fseqlen, length(seqs))
    #Z = Array{fseqlen, length(seqs)}
    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1
    close(f)
    return Array{Int8,2}(Z')
end

let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
    # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end

