module ToolsDCA


##############################
#some tools for DCA analysis
##############################
#N = num of spins (proteins)
#q = number of colors (residues, q=20+1  for proteins)
#N.b. the input parameters for gauge transformation and energy calculation must be of dimension h=(q,N) and J=(q,q,N,N)

export compute_energy_single_sequence,
    compute_energy_MSA, 
    compute_single_mut,
    fasta2matrix, 
    ind_shuffle,
    switch_gauge,
    compute_PPV,
    plot_PPV,
    compute_ROC,
    plot_ROC


include("compute_energy.jl")
include("mut_landscape.jl")
include("fasta_to_matrix.jl")
include("switch_gauge.jl")
include("PPV_ROC.jl")
end

