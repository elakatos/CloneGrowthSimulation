using Distributions, StatsBase, DataFrames, GLM

type cancercell
    mutations::Array{Int64,1}
    fitness::Int64
end

function newmutations(cancercell, mutID, p)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1
    
    neoep = rand()<p
    if neoep
        cancercell.fitness = 0
    end

    return cancercell, mutID, neoep
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness))
end

function start_population(p, initial_mut, allowNeoep=true)
    mutID = 1
    N = 1
    cells = cancercell[]
    neoep_muts = Int64[]
    push!(cells,cancercell([],1))
    for i=1:initial_mut
        cells[1],mutID,neoep = newmutations(cells[1],mutID, p)
        if neoep
            push!(neoep_muts, mutID-1)
        end
    end

    if allowNeoep==false #overwrite to make sure first cell contains no neoepitopes
        cells[1].fitness=1
        neoep_muts = Int64[]
    end

    nonimm = cells[1].fitness

    return cells, mutID, neoep_muts, nonimm, N

end

function birthdeath_neoep(b0, d0, b_, d_, Nmax, p, initial_mut=10, mu=1)

    Rmax = b0+d_ #Rmax is given by d_ as it is always >= d0

    #initialize arrays and parameters
    cells, mutID, neoep_muts, nonimm, N = start_population(p, initial_mut, false)
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
        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell
        Nt = N
        
        if cells[randcell].fitness==0 #set death rate according to whether cell is antigenic or not
            d = d_
        else
            d = d0
        end

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - cells[randcell].fitness

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep = newmutations(cells[randcell],mutID, p)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep = newmutations(cells[end],mutID, p)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + cells[randcell].fitness + cells[end].fitness
            
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
            nonimm = nonimm - cells[randcell].fitness

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
            cells, mutID, neoep_muts, nonimm, N = start_population(p, initial_mut, false)
            push!(Nvec,N)
            nonimm = cells[1].fitness
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, neoep_muts, cells, nonimmvec
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


for i=1:1000
    Nvec, tvec, mutID, neoep_muts, cells, immune = birthdeath_neoep(1, d0, 1, d_, popSize, p, initial_mut, mu);
    outNDFsim = DataFrame(t=tvec[end], N=Nvec[end], nonImm=immune[end])
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size at the end of simulation

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

    writedlm("neoep_mutations_"*string(i)*".txt", neoep_muts) #Save which mutations are neoepitope ones

end

