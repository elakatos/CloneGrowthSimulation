#=
    Julia code implementing simulations of a stochastic branching process model of tumour growth and mutation accumulations, under negative selection from the immune system due to randomly arising antigenic mutations. This code corresponds to an ideal population/measurement, with no measurement noise (but intrinsic stochasticity) and user-defined detection limit for mutations.

    For details of the model see Lakatos et al, biorXiv, 2019: https://www.biorxiv.org/content/10.1101/536433v1 and the Readme at https://github.com/elakatos/CloneGrowthSimulation

    This file implements the base simulations with fixed immune escape introduced at a preset population size and escaped and non-escaped subclones processed independently.
    Parameters: d0, popSize, detLim, p, initial_mut, mu, neoep_dist, pesc
    Outputs:
        - all_mutations_<i>.txt / all_mutations_IT_<i>.txt / all_mutations_IT2_<i>.txt: every mutation and its antigenicity value in two tab-separated columns, for mutations before escape / escaped subclone mutations / non-escaped subclone mutations
        - vaf_preIT_<i>.txt / vaf_postIT_<i>.txt / vaf_postIT2_<i>.txt : every mutation present in >detLim cells and the number of cells its present, in two tab-separated columns, for mutations before escape / escaped subclone mutations / non-escaped subclone mutations

    @Author: Eszter Lakatos (e.lakatos@qmul.ac.uk)
    Inspired by Marc J. Williams' CancerSeqSim (https://github.com/marcjwilliams1/CancerSeqSim.jl)
=#

using Distributions, StatsBase, DataFrames, GLM

function getFitness(n)
        (1 + s*n)
end

type cancercell
    mutations::Array{Int64,1}
    fitness::Float64
    epnumber::Float64
    escaped::Bool
end

function newmutations(cancercell, mutID, p, neoep_dist)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1
    
    neoep = rand()<p
    if neoep
        neoep_value = max.(0, rand(neoep_dist) ) #Take distribution from top-defined params
        cancercell.epnumber = cancercell.epnumber + neoep_value
        cancercell.fitness = getFitness(cancercell.epnumber) #fitness is affected by the number of mutations
    else
        neoep_value = 0
    end

    return cancercell, mutID, neoep_value
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber), copy(cancercellold.escaped))
end

function start_population(p, neoep_dist, initial_mut, allowNeoep=true, immThresh=0.5)
    mutID = 1
    N = 1
    cells = cancercell[]
    muts = Dict{Int64, Float64}()
    push!(cells,cancercell([],1,0,false))
    for i=1:initial_mut
        if allowNeoep
            cells[1],mutID,neoep_val = newmutations(cells[1],mutID, p, neoep_dist)
        else
            cells[1],mutID,neoep_val = newmutations(cells[1],mutID, 0, neoep_dist)
        end

        muts[mutID-1] = neoep_val
    end

    nonimm = 1*(cells[1].epnumber<immThresh)

    return cells, mutID, muts, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p, neoep_dist, initial_mut=10, mu=1, immThresh=0.5)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, muts, nonimm, N = start_population(p, neoep_dist, initial_mut)
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimmvec = Int64[]
    push!(nonimmvec, nonimm)

    while (N < Nmax) & (t < 300) #set so we can exit simulation where there is a lot of death

        #pick a random cell
        randcell = rand(1:N)
        Nt = N
        
        #a cell's immunogenicity depends on its fitness, i.e. the summed antigenicity of neoepitopes
        d = max(0, (d0 - b0)*cells[randcell].fitness + b0)

        if cells[randcell].escaped # discard effect of fitness if the cell escape
            d = d0
        end

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell     

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(cells[randcell].epnumber<immThresh)

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep_val = newmutations(cells[randcell],mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep_val = newmutations(cells[end],mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(cells[randcell].epnumber<immThresh) + 1*(cells[end].epnumber<immThresh)
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        #if r has neither birth or death (only possible if it is a non-immunogenic cell), nothing happens
        if  (b0+d)<= r
          push!(Nvec, N)
          push!(nonimmvec,nonimm)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if r > b but < d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - 1*(cells[randcell].epnumber<immThresh)

            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #if every cell dies, restart simulation from a single cell again
        if (N == 0)
            cells, mutID, muts, nonimm, N = start_population(p, neoep_dist, initial_mut)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, muts, cells, nonimmvec
end

function tumourgrow_post_it(b0, d0, Tmax, p, neoep_dist, cells, mutID, muts, mu=1, immThresh=0.5)
    #Function to track the population after some change in its environment or structure

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters: use the cells inherited from prior simulation to define quantities
    startPopSize = length(cells)
    N = startPopSize
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimm = 0
    for i=1:length(cells)
        nonimm = nonimm + 1*(cells[i].epnumber<immThresh)
    end
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while t < Tmax

        #pick a random cell
        randcell = rand(1:N)
        Nt = N
        
        #a cell's immunogenicity depends on its fitness, i.e. the summed antigenicity of neoepitopes
        d = max(0, (d0 - b0)*cells[randcell].fitness + b0)

        if cells[randcell].escaped # discard effect of fitness if the cell escape
            d = d0
        end

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell     

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(cells[randcell].epnumber<immThresh)

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep_val = newmutations(cells[randcell],mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep_val = newmutations(cells[end],mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(cells[randcell].epnumber<immThresh) + 1*(cells[end].epnumber<immThresh)
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        #if r has neither birth or death (only possible if it is a non-immunogenic cell), nothing happens
        if  (b0+d)<= r
          push!(Nvec, N)
          push!(nonimmvec,nonimm)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if r > b but < d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - 1*(cells[randcell].epnumber<immThresh)

            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies, automatically jump to Tmax
        if (N == 0)
          t = Tmax
        end

        #if we reach a high population showing exponential growth, also jump to Tmax
        if (N >= 50000)
         t = Tmax
        end

    end
    
    return Nvec, tvec, cells,nonimmvec, muts, mutID
end

function process_mutations(cells, detLim)
    mutVec = []
    for i=1:length(cells)
        append!(mutVec, cells[i].mutations)
    end

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))

    println("Mutations processed for ", length(cells), " cells.")
    return detMutDict

end


for i=1:50
    Nvec, tvec, mutID, muts, cells, immune = birthdeath_neoep(1, d0, popSize, p, neoep_dist, initial_mut, mu);

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

    writedlm("all_mutations_"*string(i)*".txt", muts) #Output dictionary storing mutations and their immunogenicity

    # Go through cell population and divide it up to two subclones, which gets grown afterwards
    subcloneSize = convert(Int64, floor(popSize/2))
    popSize = convert(Int64, popSize)

    for j=1:subcloneSize
        cells[j].escaped = true
    end

    cells_sc1 = cells[1:subcloneSize] # escaped subclone
    cells_sc2 = cells[(subcloneSize+1):popSize] # not escaped subclone
    muts_clonal = copy(muts)

    # Run simulation for both subclones from this time
    NvecIT, tvecIT, cellsIT, immuneIT, mutsIT, mutIDIT = tumourgrow_post_it(1, d0, Tmax, p, neoep_dist, cells_sc1, mutID, muts, mu)
    detMutDictIT = process_mutations(cellsIT, detLim)
    writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT)
    writedlm("all_mutationsIT_"*string(i)*".txt", mutsIT)

    NvecIT2, tvecIT2, cellsIT2, immuneIT2, mutsIT2, mutIDIT2 = tumourgrow_post_it(1, d0, Tmax, p, neoep_dist, cells_sc2, mutIDIT, muts_clonal, mu)
    detMutDictIT2 = process_mutations(cellsIT2, detLim)
    writedlm("vaf_postIT2_"*string(i)*".txt",detMutDictIT2)
    writedlm("all_mutationsIT2_"*string(i)*".txt", mutsIT2)

end

