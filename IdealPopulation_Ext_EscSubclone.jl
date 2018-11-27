using Distributions, StatsBase, DataFrames, GLM

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
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber), copy(cancercellold.escaped))
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
        neoep_muts = Int64[]
    end

    nonimm = 1*(cells[1].fitness==1)

    return cells, mutID, neoep_muts, esc_muts, nonimm, N

end

function birthdeath_neoep(b0, d0, b_, d_, Nmax, p, initial_mut, mu, pesc)

    Rmax = b0+d_ #Rmax is given by d_ as it is always >= d0

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
        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell
        Nt = N
        
        if cells[randcell].escaped # discard effect of fitness if the cell escape
                d = d0
        else
                d = d0 + (1-cells[randcell].fitness) * (d_ - d0) #otherwise set death rate accordingly
        end

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(cells[randcell].fitness==1)

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
            nonimm = nonimm + 1*(cells[randcell].fitness==1) + 1*(cells[end].fitness==1)
            
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
            nonimm = nonimm - 1*(cells[randcell].fitness==1)

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


function tumourgrow_post_it(b0, d0, b_, d_, Tmax, cells, mutID, neoep_muts, esc_muts, mu, pesc)
    #Function to track the population after it was subject to (immuno-)therapy
    #Follow through which mutations are neo-epitopes from previous simulation
    #d_ and its prevalence is modified from previous simulation according to type of therapy

    Rmax = b0+d_

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
        nonimm = nonimm + 1*(cells[i].fitness==1)
    end
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while t < Tmax #run simulation until a set time instead of population in this case

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell
        Nt = N
        
        if cells[randcell].escaped # discard effect of fitness if the cell escape
                d = d0
        else
                d = d0 + (1-cells[randcell].fitness) * (d_ - d0) #otherwise set death rate accordingly
        end

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1*(cells[randcell].fitness==1)

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
            nonimm = nonimm + 1*(cells[randcell].fitness==1) + 1*(cells[end].fitness==1)
            
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
            nonimm = nonimm - 1*(cells[randcell].fitness==1)

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
       if (N > 1.5*startPopSize)
         t = Tmax
       end

    end
    
    return Nvec, tvec, cells,nonimmvec, neoep_muts, esc_muts
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

simNum = 0
for i=1:50
    Nvec, tvec, mutID, neoep_muts, cells, immune, esc_muts = birthdeath_neoep(1, d0, 1, d_, popSize, p, initial_mut, mu, pesc);

    detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs
    writedlm("neoep_mutations_"*string(i)*".txt", neoep_muts) #Save which mutations are neoepitope ones

    # Go through cell population and divide it up to two subclones, which gets grown afterwards
    subcloneSize = convert(Int64, floor(popSize/2))
    popSize = convert(Int64, popSize)

    for j=1:subcloneSize
        cells[j].escaped = true
    end

    cells_sc1 = cells[1:subcloneSize] # escaped subclone
    cells_sc2 = cells[(subcloneSize+1):popSize] # not escaped subclone
    neoep_muts_clonal = copy(neoep_muts)

    # Run simulation for both subclones from this time
    NvecIT, tvecIT, cellsIT, immuneIT, neoep_mutsIT, esc_mutsIT = tumourgrow_post_it(1, d0, 1, d_, Tmax, cells_sc1, mutID, neoep_muts, esc_muts, mu, pesc)
    detMutDictIT = process_mutations(cellsIT, detLim)
    writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT)

    writedlm("neoep_mutationsIT_"*string(i)*".txt", neoep_mutsIT)

    NvecIT2, tvecIT2, cellsIT2, immuneIT2, neoep_mutsIT2, esc_mutsIT2 = tumourgrow_post_it(1, d0, 1, d_, Tmax, cells_sc2, mutID, neoep_muts_clonal, esc_muts, mu, pesc)
    detMutDictIT2 = process_mutations(cellsIT2, detLim)
    writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT2)

    writedlm("neoep_mutationsIT2_"*string(i)*".txt", neoep_mutsIT2)

    end

end

