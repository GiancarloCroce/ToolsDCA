function switch_gauge(h0::Array{Float64,2},
                      J0::Array{Float64,4},
                      mode::AbstractString)


    if mode=="latticeGas"
        println("Switching to lattice gas gauge")

        q, N = size(h0)

        println("N=",N,", q=",q) 

        J = zeros(J0)
        h = zeros(h0)

        #switch gauge for couplings
        for i = 1:N
            for j = 1:N
                J[:,:,i,j] = J0[:,:,i,j] .- repmat(J0[:,q,i,j],1,q) .- repmat(J0[q,:,i,j]',q,1) .+  J0[q,q,i,j] 
            end
        end

        #switch gauge for fields
        for i = 1:N
            h[:,i] = h0[:,i] .- h0[q,i]
            for j = 1:N
                if( j != i)
                    h[:,i] = h[:,i] .+ (J0[:,q,i,j] .- J0[q,q,i,j])
                end
            end
        end

        ##do element by elemnt, slow, but ok to check -> everything OK

        #Jt = zeros(J0)
        #ht = zeros(h0)

        ##switch gauge for couplings
        #for i = 1:N
        #    for j = 1:N
        #        for a = 1:q
        #            for b = 1:q
        #                Jt[a,b,i,j] = J0[a,b,i,j] -  J0[a,q,i,j] - J0[q,b,i,j] + J0[q,q,i,j]
        #            end
        #        end

        #    end
        #end

        ##switch gauge for fields
        #for i = 1:N
        #    for a = 1:q
        #        ht[a,i] = h0[a,i] - h0[q,i] 
        #        for j = 1:N
        #            if( j != i)
        #                ht[a,i] += J0[a,q,i,j] - J0[q,q,i,j]
        #            end
        #        end
        #    end
        #end

        return h,J#,ht,Jt




    elseif mode=="0sum"
        println("Switching to zero sum gauge, maximize Frobieus norm")

        q, N = size(h0)

        println("N=",N,", q=",q) 

        J = zeros(J0)
        h = zeros(h0)

        #switch gauge for couplings
        for i = 1:N
            for j = 1:N
                J[:,:,i,j] = J0[:,:,i,j].-repmat(mean(J0[:,:,i,j],1),q,1) - repmat(mean(J0[:,:,i,j],2),1,q) .+ mean(J0[:,:,i,j])
            end
        end

        #switch gauge for fields
        h[:,:] = h0[:,:] .- repmat(mean(h0,1),q) 
        for i = 1:N
            for j = 1:N
                if( j != i)
                    h[:,i] = h[:,i] .+ ( squeeze(mean(J0[:,:,i,j],2),2) .- mean(J0[:,:,i,j]) )
                end
            end
        end



        ##do element by elemnt, slow, but ok to check -> everything OK
        ##TRY!!!!

        #Jt = zeros(J0)
        #ht = zeros(h0)

        ##switch gauge for couplings
        #for i = 1:N
        #    for j = 1:N
        #        for a = 1:q
        #            for b = 1:q
        #                Jt[a,b,i,j] = J0[a,b,i,j] -  mean(J0[a,:,i,j]) - mean(J0[:,b,i,j]) + mean(J0[:,:,i,j])
        #            end
        #        end

        #    end
        #end

        ##switch gauge for fields
        #for i = 1:N
        #    for a = 1:q
        #        ht[a,i] = h0[a,i] - mean(h0[:,i]) 
        #        for j = 1:N
        #            if( j != i)
        #                ht[a,i] += mean(J0[a,:,i,j]) - mean(J0[:,:,i,j])
        #            end
        #        end
        #    end
        #end

        return h,J#,ht,Jt


    elseif (mode!="0sum" || mode !="latticeGas")
        error("Invalid mode, only accepted: 0sum or latticeGas")
    end



end



#=

################################################################################################################################
##N.B. The lattice gas switch is the only one that i checked the others probably do not work.... (problems going from MATLAB to Julia with repmat and transpose)
    if mode == "LatticeGas"
        println("Switching to lattice gas gauge")
        N = Int(maximum(size(hi))/q)
        ho = zeros(1,q*N)
        Jo = zeros(q*N,q*N)

        Jo = Ji
        ho = hi
        for i = 1:N
            for j = 1:N
                Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] - repmat(Ji[(i-1)*q+(1:q),(j-1)*q+q],1,q)
                Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] - repmat(Ji[(i-1)*q+q,(j-1)*q+(1:q)]',q,1)
                Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] =
                Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] + ones(q,q)*Ji[(i-1)*q+q,(j-1)*q+q]
                ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] + Ji[(i-1)*q+(1:q),(j-1)*q+q]
                ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] - ones(q)*(Ji[(i-1)*q+q,(j-1)*q+q])
            end
        end

        for i = 1:N
            ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] - ones(q)*(hi[(i-1)*q+q]) #repmat(ho[(i-1)*q+q],1,q)
        end
        # ho[1]=1000
        ################################################################################################################################

    elseif mode == "regl2"
        println("Switching to lattice regl2 gauge (gauge that we get in plmDCA with l2 regularization..see paper)")
        N = Int(maximum(size(hi))/q)
        ho = zeros(1,q*N)
        Jo = zeros(q*N,q*N)

        hi,Ji = switch_gauge(Ji,hi,"0sum",q)
        u = hi / (q+N-1)

        for i = 1:N
            ho[(i-1)*q+(1:q)] = hi[(i-1)*q+(1:q)] - (N-1)*u[(i-1)*q+(1:q)]
            for j = (i+1):N
                Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)] + repmat(u[(i-1)*q+(1:q)]',1,q) + repmat(u[(j-1)*q+(1:q)],q,1)
            end
        end
        Jo = Jo + Jo'


    elseif mode == "0sum"
        println("Switching to zero sum gauge, maximize Frobieus norm")

        # Works for q = 2 at l
            #       east -- and q = 21

            N = Int(maximum(size(hi))/q)
            ho = zeros(1,q*N)
            Jo = zeros(q*N,q*N)

            for i = 1:N
                for j = (i+1):N
                    Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)] - repmat(mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)],1),q,1) - repmat(mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)],2),1,q) + mean(mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)]))
                    # ho[(i-1)*q+(1:q)] =
                    # hi[(i-1)*q+(1:q)] +
                    # mean(Ji
                    # [(i-1)*q+(1:q),(j-1)*q+(1:q)],2)'
                    # -
                    # mean(mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)]))
                    # ho[(j-1)*q+(1:q)] =
                    # hi[(i-1)*q+(1:q)] +
                    # mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)],1)
                end
            end
            Jo = Jo + Jo'

            for i = 1:N
                ho[(i-1)*q+(1:q)] = hi[(i-1)*q+(1:q)]
                for j = 1:N
                    ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] + mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)],2)' - mean(mean(Ji[(i-1)*q+(1:q),(j-1)*q+(1:q)]))
                end
                ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] - mean(hi[(i-1)*q+(1:q)])
            end

        elseif mode == "wt"
            N = Int(maximum(size(hi))/q)
            Jo = Ji
            ho = hi
            for i = 1:N
                alpha = find(wt01[(i-1)*q+(1:q)]==1)
                for j = 1:N
                    beta = find(wt01[(j-1)*q+(1:q)]==1)
                    Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] - repmat(Ji[(i-1)*q+(1:q),(j-1)*q+beta()],1,q)
                    Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] - repmat(Ji[(i-1)*q+alpha(),(j-1)*q+(1:q)],q,1)
                    Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] = Jo[(i-1)*q+(1:q),(j-1)*q+(1:q)] + repmat(Ji[(i-1)*q+alpha(),(j-1)*q+beta()],q,q);
                    ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] + Ji[(i-1)*q+(1:q),(j-1)*q+beta()]'
                    ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] - repmat(Ji[(i-1)*q+alpha(),(j-1)*q+beta()],1,q)
                end
            end

            for i = 1:N
                alpha = find(wt01[(i-1)*q+(1:q)]==1)
                ho[(i-1)*q+(1:q)] = ho[(i-1)*q+(1:q)] - repmat(hi[(i-1)*q+alpha()],1,q) #repmat(ho[(i-1)*q+q],1,q)
            end
        end

        return ho,Jo
    end

    =#
