#function to compute the energy of the single-double mutants 

function compute_single_mut(h::Array{Float64,2},
                            J::Array{Float64,4},
                            wild_type::AbstractString)

alphabet = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21]


    wt = fasta2matrix(wild_type,0.9)
    e0 = compute_energy_single_sequence(h,J,wt[:])
    _,N = size(wt)

    res = open("sm_"*wild_type, "w")

    for k = 1:N
        for i in alphabet
            tmp = copy(wt)
            tmp[k] = i
            e_mut = compute_energy_single_sequence(h,J,tmp[:])
            println(res,k," ",i," ", e_mut, " ", e_mut - e0)
        end
    end
    close(res)
end


