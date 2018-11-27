using Distributions, StatsBase, DataFrames, GLM

type cancercell
    mutations::Array{Int64,1}
    fitness::Float64
    epnumber::Int64
    escaped1::Bool #Two types of immune escape documented separately
    escaped2::Bool
end

function newmutations(cancercell, mutID, p, pesc, esc_effect)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1

    mut_escape1 = rand()<pesc[1]
    mut_escape2 = rand()<pesc[2]
    if mut_escape1
        cancercell.escaped1 = true
    end
    if mut_escape2
        cancercell.escaped2 = true
        cancercell.fitness = getFitness(convert(Int64,floor(cancercell.epnumber * (1-esc_effect) )))
    end

    neoep = rand()<p
    if neoep
        cancercell.epnumber = cancercell.epnumber + 1
        if cancercell.escaped2
            cancercell.fitness = getFitness(convert(Int64,floor(cancercell.epnumber * (1-esc_effect) )))
        else
            cancercell.fitness = getFitness(cancercell.epnumber) #fitness is affected by the number of mutations
        end
    end

    return cancercell, mutID, neoep, mut_escape1, mut_escape2
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber), copy(cancercellold.escaped1), copy(cancercellold.escaped2))
end

function start_population(p, initial_mut, pesc, esc_effect, allowNeoep=true)
    mutID = 1
    N = 1
    cells = cancercell[]
    neoep_muts = Int64[]
    esc_muts = Dict{Int64, Int64}()
    push!(cells,cancercell([],1, 0, false, false))
    for i=1:initial_mut
        cells[1],mutID,neoep, mut_escape1, mut_escape2 = newmutations(cells[1],mutID, p,pesc, esc_effect)
        if neoep
            push!(neoep_muts, mutID-1)
        end
        if mut_escape1
            esc_muts[mutID-1]=1
        end
        if mut_escape2
            esc_muts[mutID-1]=2
        end
    end

    if allowNeoep==false #overwrite to make sure first cell contains no neoepitopes
        cells[1].fitness=1
        neoep_muts = Int64[]
    end

    # Compute non-immune as: 1 for any without neoantigens. esc_effect for any with neoantigen, but with immune escape.
    # 0 for first type of immune escape (PD-L1-type)
    nonimm = 1*(cells[1].fitness==1) + esc_effect*(cells[1].fitness<1 && cells[1].escaped2==true)

    return cells, mutID, neoep_muts, esc_muts, nonimm, N

end

function birthdeath_neoep(b0, d0, b_, d_, Nmax, p, initial_mut, mu, pesc, esc_effect)

    Rmax = b0+d_ #Rmax is given by d_ as it is always >= d0

    #initialize arrays and parameters
    cells, mutID, neoep_muts, esc_muts, nonimm, N = start_population(p, initial_mut, pesc, esc_effect)
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while (N < Nmax) & (t < 300) #set so we can exit simulation where there is a lot of death

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell
        Nt = N
        
        if cells[randcell].escaped1 #discard effect of fitness if escaped
                d = d0
        else
            #Cells with type2 escape and without should be computed according to fitness that takes neoep number and esc into account
            d = d0 + (1-cells[randcell].fitness) * (d_ - d0)
        end

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - (1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true))

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep,mut_escape1, mut_escape2 = newmutations(cells[randcell],mutID, p, pesc, esc_effect)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape1
                    esc_muts[mutID-1]=1
                    print("Escape type 1 at ", t, " time.")
                end
                if mut_escape2
                    esc_muts[mutID-1]=2
                    print("Escape type 2 at ", t, " time.")
                end
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep,mut_escape1, mut_escape2 = newmutations(cells[end],mutID, p, pesc, esc_effect)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape1
                    esc_muts[mutID-1]=1
                end
                if mut_escape2
                    esc_muts[mutID-1]=2
                end
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true) + 1*(cells[end].fitness==1) + esc_effect*(cells[end].fitness<1 && cells[end].escaped2==true)
            
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
            nonimm = nonimm - (1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true))

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
            cells, mutID, neoep_muts, esc_muts, nonimm, N = start_population(p, initial_mut, pesc, esc_effect)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, neoep_muts, cells, nonimmvec, esc_muts
end

function tumourgrow_post_it(b0, d0, b_, d_, Tmax, cells, mutID, neoep_muts, esc_muts, mu, pesc, esc_effect)
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
        nonimm = nonimm + 1*(cells[i].fitness==1) + esc_effect*(cells[i].fitness<1 && cells[i].escaped2==true)
    end
    nonimmvec = Float64[]
    push!(nonimmvec, nonimm)

    while t < Tmax #run simulation until a set time instead of population in this case

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        #set death rate according to fitness, escape status does not matter due to nature of therapy
        d = d0 + (1-cells[randcell].fitness) * (d_ - d0)

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - (1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true))

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            #still record escape mutations, but type1 does not carry any advantage anymore
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID,neoep,mut_escape1, mut_escape2 = newmutations(cells[randcell],mutID, p, pesc, esc_effect)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape1
                    esc_muts[mutID-1]=1
                    print("Escape type 1 at ", t, " time.")
                end
                if mut_escape2
                    esc_muts[mutID-1]=2
                    print("Escape type 2 at ", t, " time.")
                end
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID,neoep,mut_escape1, mut_escape2 = newmutations(cells[end],mutID, p, pesc, esc_effect)
                if neoep
                    push!(neoep_muts, mutID-1)
                end
                if mut_escape1
                    esc_muts[mutID-1]=1
                end
                if mut_escape2
                    esc_muts[mutID-1]=2
                end
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true) + 1*(cells[end].fitness==1) + esc_effect*(cells[end].fitness<1 && cells[end].escaped2==true)
            
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
            nonimm = nonimm - (1*(cells[randcell].fitness==1) + esc_effect*(cells[randcell].fitness<1 && cells[randcell].escaped2==true))

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
    escPopulations = Dict{Int64, Int64}()
    escPopulations[1] = 0; escPopulations[2] = 0; escPopulations[3] = 0;
    for i=1:length(cells)
        append!(mutVec, cells[i].mutations)
        if cells[i].escaped1
            if cells[i].escaped2
                escPopulations[3] = escPopulations[3]+1
            else
                escPopulations[1] = escPopulations[1]+1
            end
        else
            if cells[i].escaped2
                escPopulations[2] = escPopulations[2]+1
            end
        end

    end

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))

    println("Mutations processed for ", length(cells), " cells.")
    return detMutDict, escPopulations

end

simNum = 0
while (simNum < 100)
    Nvec, tvec, mutID, neoep_muts, cells, immune, esc_muts = birthdeath_neoep(1, d0, 1, d_, popSize, p, initial_mut, mu, pesc, esc_effect);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)

    detMutDict, escPopulations = process_mutations(cells, detLim)
    #Only proceed if there are any escaped cells
    if (escPopulations[1]+escPopulations[2]+escPopulations[3])>0
        simNum = simNum + 1
        i = simNum
        writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation, but only in escaped cases
        writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs

        writedlm("neoep_mutations_"*string(i)*".txt", neoep_muts) #Save which mutations are neoepitope ones
        writedlm("escape_mutations_"*string(i)*".txt", esc_muts) #Write mutations (if they caused escape) and type of escape
        writedlm("escaped_populations_"*string(i)*".txt", escPopulations) #Write number of cells with each escape type)

        if doTherapy
            #Apply immunotherapy to the population and record information afterwards
            NvecIT, tvecIT, cellsIT, immuneIT, neoep_mutsIT, esc_mutsIT = tumourgrow_post_it(1, d0, 1, d_it, Tmax, cells, mutID, neoep_muts, esc_muts, mu, pesc, esc_effect)
            outNDFsimIT = DataFrame(t=tvecIT, N=NvecIT, nonImm=immuneIT)
            detMutDictIT, escPopulationsIT = process_mutations(cellsIT, detLim)
            writetable("postIT_"*string(i)*".txt", outNDFsimIT)
            writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT)

            writedlm("neoep_mutationsIT_"*string(i)*".txt", neoep_mutsIT)
            writedlm("escape_mutationsIT_"*string(i)*".txt", esc_mutsIT)
            writedlm("escaped_populationsIT_"*string(i)*".txt", escPopulationsIT)
        end

    end

end

