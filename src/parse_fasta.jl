using FastaIO
using Compat

export parse_fasta 

function parse_fasta(filename::AbstractString, max_gap_fraction::Real)

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

    parsed_fasta = open("reduced_"*filename, "w")

    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue

        parsed_seq = ""#Array{Any}(fseqlen)
        for i = 1:fseqlen
            c = seq[inds[i]]

            #write reduced fasta
            parsed_seq *= string(c)
        end
        seqid += 1
        @printf(parsed_fasta, "%s\n", name )
        @printf(parsed_fasta, "%s\n", parsed_seq)
    end
    @assert seqid == length(seqs) + 1

    close(f)
    close(parsed_fasta)

end
