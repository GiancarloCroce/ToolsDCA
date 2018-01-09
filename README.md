# some tools for DCA

## install: 
    julia> Pkg.clone("https://github.com/GiancarloCroce/ToolsDCA")

## convections
* N = num of spins
* M = num of sequences(sampled configurations)
* q = num of colors for (generalized) Potts model
* fields must have dimension: h=(q,N)
* couplings must have dimension: J=(q,q,N,N)

### list of implemented functions:

   - **compute_energy.jl**: Energy calculation (for mutational landscape studies):
        - compute_energy_single_sequence(h,J,S::vector of configuration)
        - compute_energy_MSA(h,J,S::matrix MxN [matrix form of a MSA])
            the mapping fasta_MSA -> matrix is the same of GaussDCA
        - compute_energy_MSA(h,J,fastafile:: a fasta file)
    
    - **fasta_to_matrix.jl**: map MSA in matrix and shuffling for indipendent model
        - fasta2matrix: map fasta MSA in matrix (same of GaussDCA, see "alphabet")
        - shuffle_matrix(matrix, dim, num_of_shuffle) 
            matrix: the matrix to be shuffled
            if( dim == 0): shuffle all elements, no statistics is conserved
            if( dim == 1): for each column, indipendently, shuffle along rows (conserv 1-point frequencies, destroy 2-points correlation) [i.e. inferred J_ij will be close to 0 ( zero + finite size and phylogenetic effect) ]
            if( dim == 2): for each row, indipendently, shuffle along columns: keep the magnetization of the row

    - **switch_gauge.jl**: to switch along different gauge
        - switch_gauge(h,J,mode), mode can be "lattice gas gauge" and "0sum gauge" (max Froebius norm).
            Note: if you use plmDCA ASYMMETRIC the gauge for h_i should be changed indipendently for h_i coming from J_ix and for h_i coming from J_xi (where x is everything exept i)  ---> switch gauge is not good for (h,J) coming out of plmDCA  

    - **PPV_ROC.jl**: compute and plot PPV and ROC curve
        - compute_PPV(v::vector of tp [already sorted according to some measurem eg. Froebius norm of DCA couplings)
        - compute_ROC(v:: as above)
        - plot_ROC(v::as above): plot ROC curve and compute the AUC 
            
