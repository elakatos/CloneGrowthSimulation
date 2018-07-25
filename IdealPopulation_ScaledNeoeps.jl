using Distributions, StatsBase, DataFrames, GLM

type cancercell
    mutations::Array{Int64,1}
    fitness::Float64
end

function newmutations(cancercell, mutID, p, neoep_dist)

    numbermutations= 1
    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    
    neoep = rand()<p
    if neoep
        neoep_value = min.(1, rand(neoep_dist)/10 ) #Take distribution from top-defined params
        cancercell.fitness = min.(1-neoep_value, cancercell.fitness) #Take the minimum as fitness value
    else
        neoep_value = 0
    end

    return cancercell, mutID, neoep_value
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness))
end

function tumourgrow_birthdeath_neoep(b0, d0, b_, d_, Nmax, p, neoep_dist)


    #Rmax is defined from the highest possible death and birthrate
    Rmax = b0+d_

    #initialize arrays and parameters
    mutID = 1
    cells = cancercell[]
    muts = Dict{Int64, Float64}()
    push!(cells,cancercell([],1))
    for i=1:10
        cells[1],mutID,neoep_val = newmutations(cells[1],mutID, p, neoep_dist)
        muts[mutID] = neoep_val

    end
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimm = cells[1].fitness
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        #Set deathrate according to immunogenicity of cell preserved in fitness (fitness=1: nonimm, 0:very imm)
        d = d0 + (d_-d0)*(1-cells[randcell].fitness)

        #birth event if r<birthrate (birthrate is fixed for all cell types)
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - cells[randcell].fitness
            #add new mutations to both new cells
            cells[randcell],mutID,neoep_val = newmutations(cells[randcell],mutID, p, neoep_dist)
            muts[mutID] = neoep_val

            cells[end],mutID,neoep = newmutations(cells[end],mutID, p, neoep_dist)
            muts[mutID] = neoep_val

            #note down non/immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + cells[randcell].fitness + cells[end].fitness
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        if  (b0+d)<= r
          push!(Nvec, N)
          push!(nonimmvec,nonimm)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - cells[randcell].fitness
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies reinitialize simulation
        if (N == 0)
                mutID = 1
        cells = cancercell[]
        muts = Dict{Int64, Float64}()
        push!(cells,cancercell([],1))
        for i=1:10
            cells[1],mutID,neoep_val = newmutations(cells[1],mutID, p, neoep_dist)
            muts[mutID] = neoep_val

        end
        N = 1
        push!(Nvec,N)
            nonimm = cells[1].fitness
        push!(nonimmvec, nonimm)
        push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, muts, cells, nonimmvec
end

function tumourgrow_post_it(b0, d0, b_, d_, Tmax, startPopSize, cells, mutID, muts, p, neoep_dist)
    #Function only to track the size of population afterwards, don't care about mutants
    #Either immunogenicity status is fixed and heritable or new mutations can change it
    #Follow through which mutations are neo-epitopes from previous simulation

    #Rmax starts with b + d
    Rmax = b0+d_

    #initialize arrays and parameters: cells are inputed, only time and pop-size matters
    startPopSize = length(cells)
    N = startPopSize
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimm = 0
    for i=1:length(cells)
        nonimm = nonimm + cells[i].fitness
    end
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while t < Tmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        #Set deathrate according to immunogenicity of cell preserved in fitness (fitness=1: nonimm, 0:very imm)
        d = d0 + (d_-d0)*(1-cells[randcell].fitness)

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - cells[randcell].fitness
            #add new mutations to both new cells
            cells[randcell],mutID,neoep_val = newmutations(cells[randcell],mutID, p, neoep_dist)
            muts[mutID] = neoep_val

            cells[end],mutID,neoep = newmutations(cells[end],mutID, p, neoep_dist)
            muts[mutID] = neoep_val
            #note down non/immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + cells[randcell].fitness + cells[end].fitness
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        if  (b0+d)<= r
          push!(nonimmvec, nonimm)
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - cells[randcell].fitness
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(nonimmvec, nonimm)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies, automatically jump to Tmax
        if (N == 0)
          t = Tmax
        end

       #if we reach a high population showing exponential growth, also jump to Tmax
       if (N > 1.5*startPopSize)
         t = Tmax
       end

    end
    
    return Nvec, tvec, cells,nonimmvec, muts, mutID
end

function process_mutations(cells, detLim, popSize)
    mutVec = []
    for i=1:length(cells)
        append!(mutVec, cells[i].mutations)
    end

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))
    VAF = collect(values(detMutDict))/(2*popSize)
    VAFtotal = Array{Float64}(0)
    append!(VAFtotal, VAF)
    println("Mutations processed for ", length(cells), "cells.")
    return detMutDict, VAFtotal

end

function analyse_vaf(VAFtotal, steps, fmin, fmax)
    steps = fmax:-0.0001:fmin
    cumcount = Array{Int64}(0)
    cumfm = Array{Float64}(0)
    invv = Array{Float64}(0)
    for i in steps
        push!(cumcount, sum(VAFtotal .>= i))
        push!(cumfm, sum(VAFtotal[VAFtotal .>= i]))
        push!(invv, 1/i-1/fmax)
    end
    cumcount = cumcount-cumcount[1]
    cumfm = cumfm - cumfm[1]
    vafDF = DataFrame(invf = invv, cumcount = cumcount)
    lmfit = lm(@formula(cumcount ~ invf+0), vafDF)
    rsq = r2(lmfit)

    return cumcount, cumfm, rsq
end



outRvec = []
outvafDF = DataFrame(invf = (1./steps - 1/fmax))
outfimDF = DataFrame(f=steps)
for i=1:200
    Nvec, tvec, mutID, muts, cells, immune = tumourgrow_birthdeath_neoep(1, d0, 1, d_, popSize, p, neoep_dist);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

    detMutDict, VAFtotal = process_mutations(cells, detLim, popSize)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs
    #Process mutations, i.e. cumulative distributions
    cumcount, cumfm, rsq = analyse_vaf(VAFtotal, steps, fmin, fmax)
    push!(outRvec, rsq)
    outvafDF[Symbol("sim_$(i)")] = cumcount
    outfimDF[Symbol("sim_$(i)")] = cumfm

    # #After all is recorded, we do 'generic' therapy
    # NvecT, tvecT, cellsT, immuneT, mutsT, mutIDT = tumourgrow_post_it(1, d_t, 1, d_t, TmaxT, popSize, cells, mutID, muts, p, neoep_dist)
    # outTDFsim = DataFrame(t=tvecT, N=NvecT, nonImm=immuneT )
    # writetable("midIT_"*string(i)*".txt", outTDFsim)
    # #Process the mutation after therapy just like before
    # popSizeT = length(cellsT) 
    # detMutDictT, VAFtotalT = process_mutations(cellsT, detLim, popSizeT)
    # writedlm("vaf_midIT_"*string(i)*".txt",detMutDictT)

    # #Then we go on to simulate the population after 'immunotherapy' and process mutations
    # NvecIT, tvecIT, cellsIT, immuneIT, mutsIT = tumourgrow_post_it(1, d0_it, 1, d_it, TmaxIT, popSizeT, cellsT, mutIDT, mutsT, p, neoep_dist)
    # outITDFsim = DataFrame(t=tvecIT, N=NvecIT, nonImm=immuneIT )
    # writetable("postIT_"*string(i)*".txt", outITDFsim)
    # popSizeIT = length(cellsIT)
    # detMutDictIT, VAFtotalIT = process_mutations(cellsIT, detLim, popSizeIT)
    # writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT) #Save mutation-VAF pairs post IT
    
    writedlm("all_mutations_"*string(i)*".txt", muts) #Dictionary storing mutations and their immunogenicity

end

writetable("Vafdf.csv",outvafDF)
writetable("Fimdf.csv",outfimDF)
writedlm("Rsqs.txt",outRvec)
