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

function newmutations(cancercell, mutID, p)
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

function start_population(p, initial_mut, clonalEscape, clonalNeoep, allowNeoep=true)
    mutID = 1
    N = 1
    cells = cancercell[]
    muts = Dict{Int64, Float64}()
    push!(cells,cancercell([],1,0,false))
    for i=1:initial_mut
        if allowNeoep
            cells[1],mutID,neoep_val = newmutations(cells[1],mutID, p)
        else
            cells[1],mutID,neoep_val = newmutations(cells[1],mutID, 0)
        end

        muts[mutID-1] = neoep_val
    end

    if clonalNeoep
        muts[mutID-1] = 0.8
        cells[1].epnumber = cells[1].epnumber + 0.8
        cells[1].fitness = getFitness(cells[1].epnumber)
    end

    cells[1].escaped = clonalEscape

    nonimm = 1*(cells[1].epnumber<immThresh)

    return cells, mutID, muts, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p, clonalEscape,clonalNeoep, initial_mut=10, mu=1)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, muts, nonimm, N = start_population(p, initial_mut, clonalEscape, clonalNeoep)
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
                cells[randcell],mutID,neoep_val = newmutations(cells[randcell],mutID, p)
                muts[mutID-1] = neoep_val
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep_val = newmutations(cells[end],mutID, p)
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
            cells, mutID, muts, nonimm, N = start_population(p, initial_mut, clonalEscape, clonalNeoep)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, muts, cells, nonimmvec
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
    Nvec, tvec, mutID, muts, cells, immune = birthdeath_neoep(1, d0, popSize, p, clonalEscape, clonalNeoep, initial_mut, mu);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

    writedlm("all_mutations_"*string(i)*".txt", muts) #Dictionary storing mutations and their immunogenicity

end

