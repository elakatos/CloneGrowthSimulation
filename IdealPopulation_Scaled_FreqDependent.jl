#=
    Julia code implementing simulations of a stochastic branching process model of tumour growth and mutation accumulations, under negative selection from the immune system due to randomly arising antigenic mutations. This code corresponds to an ideal population/measurement, with no measurement noise (but intrinsic stochasticity) and user-defined detection limit for mutations.

    For details of the model see Lakatos et al, biorXiv, 2019: https://www.biorxiv.org/content/10.1101/536433v1 and the Readme at https://github.com/elakatos/CloneGrowthSimulation

    This file extends the base simulations described above by including frequency-dependencies in immune predation (negative selection) - so that antigens present in a large number of cells have a higher weight in defining antigenicity.
    Parameters: d0, popSize, detLim, p, initial_mut, mu, neoep_dist, frequency_func
    Outputs:
        - neoep_mutations_<i>.txt : every antigenic mutations, its antigenicity and its frequency in tab-separated columns
        - preIT_<i>.txt : cell growth curve with time, total population and non-immunogenic population in three comma-separated columns
        - vaf_preIT_<i>.txt : every mutation present in >detLim cells and the number of cells its present, two tab-separated columns

    @Author: Eszter Lakatos (e.lakatos@qmul.ac.uk)
    Inspired by Marc J. Williams' CancerSeqSim (https://github.com/marcjwilliams1/CancerSeqSim.jl)
=#

using Distributions, StatsBase, DataFrames, GLM

function getFitness(cancercell, epitope_table)
    A = 0
    for i=1:length(cancercell.epitopes)
        A = A + epitope_table[cancercell.epitopes[i]][1]*frequency_func(epitope_table[cancercell.epitopes[i]][2])
    end
    return (1 + s*A)
end

type cancercell
    mutations::Array{Int64,1}
    epitopes::Array{Int64,1}
end

function newmutations(cancercell, mutID, p, epitope_table)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    
    neoep = rand()<p
    if neoep
        neoep_value = max.(0, rand(Exponential(0.2)) ) #Take distribution from top-defined params
        cancercell.epitopes = append!(cancercell.epitopes, mutID)
        epitope_table[mutID] = (neoep_value, 1)
    end
    mutID = mutID + 1

    return cancercell, mutID, epitope_table
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.epitopes))
end

function adjust_epcount(cancercell, epitope_table, adjustment)
    for i=1:length(cancercell.epitopes)
        old = epitope_table[cancercell.epitopes[i]]
        epitope_table[cancercell.epitopes[i]] = (old[1], old[2]+adjustment)
    end
    return epitope_table
end

function start_population(p, initial_mut, allowNeoep=true, immThresh=0.5)
    mutID = 1
    N = 1
    cells = cancercell[]
    epitope_table = Dict{Int64, Tuple{Float64,Int64}}()
    push!(cells,cancercell([],[]))
    for i=1:initial_mut
        if allowNeoep
            cells[1],mutID,epitope_table = newmutations(cells[1],mutID, p, epitope_table)
        else
            cells[1],mutID,epitope_table = newmutations(cells[1],mutID, 0, epitope_table)
        end
    end

    nonimm = 1*(length(cells[1].epitopes)<immThresh)

    return cells, mutID, epitope_table, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p, initial_mut=10, mu=1, immThresh=0.5)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, epitope_table, nonimm, N = start_population(p, initial_mut)
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
        # each antigen is also weighted by the number of its carriers (its frequency)
        fitness = getFitness(cells[randcell], epitope_table)
        d = max(0, (d0 - b0)*fitness + b0)

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell     

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces and increase count in all epitopes
            push!(cells, copycell(cells[randcell]))
            epitope_table = adjust_epcount(cells[randcell],epitope_table,1)

            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(length(cells[randcell].epitopes)<immThresh)

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,epitope_table = newmutations(cells[randcell],mutID, p, epitope_table)
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,epitope_table = newmutations(cells[end],mutID, p, epitope_table)
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(length(cells[randcell].epitopes)<immThresh) + 1*(length(cells[end].epitopes)<immThresh)
            
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
            nonimm = nonimm - 1*(length(cells[randcell].epitopes)<immThresh)

            #remove deleted cell and adjust epitope count
            epitope_table = adjust_epcount(cells[randcell], epitope_table, -1)
            deleteat!(cells,randcell)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #if every cell dies, restart simulation from a single cell again
        if (N == 0)
            cells, mutID, epitope_table, nonimm, N = start_population(p, initial_mut)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, epitope_table, cells, nonimmvec
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

for i=1:100
    Nvec, tvec, mutID, epitope_table, cells, immune = birthdeath_neoep(1, d0, popSize, p, initial_mut, mu);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

    writedlm("neoep_mutations_"*string(i)*".txt", epitope_table) #Output dictionary storing mutations and their antigenicity and frequency

end

