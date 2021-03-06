using Distributions, StatsBase, DataFrames, GLM

function getFitness(n)
        (1 + s*n)
end

type cancercell
    mutations::Array{Int64,1}
    fitness::Float64
    epnumber::Int64
    escaped::Bool
end

function newmutations(cancercell, mutID, p, pesc)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1
    
    neoep = rand()<p
    if neoep
        cancercell.epnumber = cancercell.epnumber + 1
        cancercell.fitness = getFitness(cancercell.epnumber) #fitness is affected by the number of mutations
    end

    mut_escape = rand()<pesc
    if mut_escape
        cancercell.escaped = true
    end

    return cancercell, mutID, neoep, mut_escape
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber),copy(cancercellold.escaped))
end

function start_population(p, initial_mut, pesc,allowNeoep=true)
    mutID = 1
    N = 1
    cells = cancercell[]
    neoep_muts = Int64[]
    esc_muts = Int64[]
    push!(cells,cancercell([],1, 0, false))
    for i=1:initial_mut
        cells[1],mutID,neoep,mut_escape = newmutations(cells[1],mutID, p,pesc)
        if neoep
            push!(neoep_muts, mutID-1)
        end
        if mut_escape
            push!(esc_muts, mutID-1)
        end
    end

    if allowNeoep==false #overwrite to make sure first cell contains no neoepitopes
        cells[1].fitness=1
        cells[1].epnumber=0
        neoep_muts = Int64[]
    end

    nonimm = 1*(cells[1].epnumber==0)

    return cells, mutID, neoep_muts, esc_muts, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p, initial_mut, mu, pesc)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, neoep_muts, esc_muts, nonimm, N = start_population(p, initial_mut, pesc)
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

        if cells[randcell].escaped # discard effect of fitness if the cell escape
                d = d0
        else
                d = (d0 - b0)*cells[randcell].fitness + b0
        end

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0
            
            #population increases by one
       if N==200
            cells[randcell].escaped = true
            println("Escape at clone size 200.\n")
        end
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(cells[randcell].epnumber==0)

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep,mut_escape = newmutations(cells[randcell],mutID, p, pesc)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape
                    push!(esc_muts, mutID-1)
                end
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep,mut_escape = newmutations(cells[end],mutID, p, pesc)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape
                    push!(esc_muts, mutID-1)
                end
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(cells[randcell].epnumber==0) + 1*(cells[end].epnumber==0)
            
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
            nonimm = nonimm - 1*(cells[randcell].epnumber==0)

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
            cells, mutID, neoep_muts, esc_muts, nonimm, N = start_population(p, initial_mut, pesc)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, neoep_muts, cells, nonimmvec, esc_muts
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


for i=1:25
    Nvec, tvec, mutID, neoep_muts, cells, immune, esc_muts = birthdeath_neoep(1, d0, popSize, p, initial_mut, mu, pesc);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

    writedlm("neoep_mutations_"*string(i)*".txt", neoep_muts) #Save which mutations are neoepitope ones
    writedlm("escape_mutations_"*string(i)*".txt", esc_muts) #Write mutations if they caused escape

end

