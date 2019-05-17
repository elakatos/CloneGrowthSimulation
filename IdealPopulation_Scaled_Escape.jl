#=
    Julia code implementing simulations of a stochastic branching process model of tumour growth and mutation accumulations, under negative selection from the immune system due to randomly arising antigenic mutations. This code corresponds to an ideal population/measurement, with no measurement noise (but intrinsic stochasticity) and user-defined detection limit for mutations.

    For details of the model see Lakatos et al, biorXiv, 2019: https://www.biorxiv.org/content/10.1101/536433v1 and the Readme at https://github.com/elakatos/CloneGrowthSimulation

    This file implements the base simulations with probabilistic active immune escape, acquired by chance throughout tumour growth.
    Parameters: d0, popSize, detLim, p, initial_mut, mu, neoep_dist, pesc
    Outputs:
        - neoep_mutations_<i>.txt : list of mutations above antigenicity threshold 0.2
        - escape_mutations_<i>.txt: list of immune escape mutations
        - preIT_<i>.txt : cell growth curve with time, total population and non-immunogenic population in three comma-separated columns
        - vaf_preIT_<i>.txt : every mutation present in >detLim cells and the number of cells its present, two tab-separated columns

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

function newmutations(cancercell, mutID, p, pesc, neoep_dist)
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

    mut_escape = rand()<pesc
    if mut_escape
        cancercell.escaped = true
    end

    return cancercell, mutID, neoep_value, mut_escape
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber),copy(cancercellold.escaped))
end

function start_population(p, initial_mut, pesc, neoep_dist, allowNeoep=true, immThresh=0.5)
    mutID = 1
    N = 1
    cells = cancercell[]
    muts = Int64[]
    esc_muts = Int64[]
    push!(cells,cancercell([],1,0, false))
    for i=1:initial_mut
        if allowNeoep
            cells[1],mutID,neoep_val,mut_escape = newmutations(cells[1],mutID, p, pesc, neoep_dist)
        else
            cells[1],mutID,neoep_val,mut_escape = newmutations(cells[1],mutID, 0, pesc, neoep_dist)
        end

        if neoep_val > 0.2
            push!(muts, mutID-1)
        end
        if mut_escape
            push!(esc_muts, mutID-1)
        end
    end

    nonimm = 1*(cells[1].epnumber<immThresh)

    return cells, mutID, muts, esc_muts, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p,pesc,neoep_dist, initial_mut=10, mu=1, immThresh=0.5)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, muts, esc_muts, nonimm, N = start_population(p, initial_mut, pesc,neoep_dist)
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
                cells[randcell],mutID,neoep_val,mut_escape = newmutations(cells[randcell],mutID, p,pesc, neoep_dist)
                if neoep_val > 0.2
                    push!(muts, mutID-1)
                end
                if mut_escape
                    push!(esc_muts, mutID-1)
                end
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep_val,mut_escape = newmutations(cells[end],mutID, p,pesc, neoep_dist)
                if neoep_val > 0.2
                    push!(muts, mutID-1)
                end
                if mut_escape
                    push!(esc_muts, mutID-1)
                end
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
            cells, mutID, muts,esc_muts, nonimm, N = start_population(p, initial_mut,pesc, neoep_dist)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, muts,esc_muts, cells, nonimmvec
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

i = 1
while (i < 51)
    Nvec, tvec, mutID, muts,esc_muts, cells, immune = birthdeath_neoep(1, d0, popSize, p,pesc,neoep_dist, initial_mut, mu);
    if Nvec[end]>0.5*popSize
        outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
        writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

        detMutDict = process_mutations(cells, detLim)
        writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

        writedlm("neoep_mutations_"*string(i)*".txt", muts) #List of neoantigen mutations
        writedlm("escape_mutations_"*string(i)*".txt", esc_muts) #Write mutations if they caused escape
        i = i+1
    end

end

