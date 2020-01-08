#=
    Julia code implementing simulations of a death-birth Moran process (fixed population size) with negative selection on immunogenic mutations.
    Parameters: maxIter, startPopSize, detLim, HLAeff (parameter of neoep distribution), p, initialAnt, mu, pesc
    Outputs:
        - time_course_<i>.txt : cell growth curve with time, total population and non-immunogenic population in three comma-separated columns
        - vaf_final_<i>.txt : every mutation present in >detLim cells and the number of cells its present, two tab-separated columns
        - neoep_mutations_<i>.txt: mutations with antigenicity value above a predefined threshold (0.2 below)

    @Author: Eszter Lakatos (e.lakatos@qmul.ac.uk)
    Inspired by Marc J. Williams' CancerSeqSim (https://github.com/marcjwilliams1/CancerSeqSim.jl)
=#

using Distributions, StatsBase, DataFrames, GLM

neoep_dist = Exponential(1/HLAeff)

function getDeath(n, s)
    b0 = 1 # b0 and d0 are the same for all cells, but taken the same as in exponentially growing populations
    d0 = 0.1
    fitness = (1 + s*n)
    return max(0, (d0 - b0)*fitness + b0)
end

type cancercell
    mutations::Array{Int64,1}
    death::Float64
    epnumber::Float64
    escaped::Bool
end

function newmutations(cancercell, mutID, p, pesc, neoep_dist, s)
    
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1
    neoep = rand()<p
    if neoep
        neoep_value = max.(0, rand(neoep_dist) ) #Take distribution from top-defined params
        cancercell.epnumber = cancercell.epnumber + neoep_value    #antigenicity value increases with mutation's antigenicity
        cancercell.death = getDeath(cancercell.epnumber, s)
    else
        neoep_value = 0
    end

    mut_escape = rand()<pesc
    if mut_escape
        cancercell.escaped = true
    end

    return cancercell, mutID, neoep_value
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.death), copy(cancercellold.epnumber),copy(cancercellold.escaped))
end


function birthdeath_neoep_constant(maxIter, popSize, p,pesc,neoep_dist, initialAnt, mu)
    
    #initialize arrays and parameters

    initialDeath = getDeath(initialAnt, s)
    cells = cancercell[]
    deathRates = Float64[]
    antValues = Float64[]
    for i = 1:popSize
        push!(cells,cancercell([],initialDeath, initialAnt, false))
        push!(deathRates, initialDeath)
        push!(antValues, initialAnt)
    end

    totAnt = Float64[]
    maxAnt = Float64[]
    sumant = initialAnt*popSize
    topant = initialAnt
    push!(totAnt, sumant)
    push!(maxAnt, topant)
    Rmax = sum(deathRates)

    muts = Int64[]
    mutID = 1
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    iter = 1

    while iter < maxIter

        iter = iter+1

        # Pick cell that dies based on death rates
        r = rand(Uniform(0,Rmax)) # random variable generated between 0 and sum of death rates
        nD = find(cumsum(deathRates).>=r)[1]
        dt = 1/(Rmax) .* - log(rand()) # time is death event and birth event

        # Pick cell that proliferates randomly from remaining cells
        nB = rand([collect(1:(nD-1)); collect((nD+1):popSize)])
        dt = dt + 1/(popSize).* - log(rand())
        t = t + dt
        
        # Replace dead cell with proliferating cell's copy and add mutations

        Rmax = Rmax - cells[nD].death - cells[nB].death # subtract old values from both summed quantities
        sumant = sumant - cells[nD].epnumber - cells[nB].epnumber
        antValues[nD] = 0

        if cells[nD].epnumber == topant
        # if selected cell was the (a) one with maximum antigenicity, re-compute maximum
            topant = maximum(antValues)
        end

        cells[nD] = copycell(cells[nB])
        for i=1:(rand(Poisson(mu)))
            cells[nD], mutID, neoep_val = newmutations(cells[nD], mutID, p,pesc, neoep_dist, s)
            if neoep_val > 0.2
                push!(muts, mutID-1)
            end
        end
        for i=1:(rand(Poisson(mu)))
            cells[nB], mutID, neoep_val = newmutations(cells[nB], mutID, p,pesc, neoep_dist, s)
            if neoep_val > 0.2
                push!(muts, mutID-1)
            end
        end

        deathRates[nD] = cells[nD].death
        deathRates[nB] = cells[nB].death
        antValues[nD] = cells[nD].epnumber
        antValues[nB] = cells[nB].epnumber
        Rmax = Rmax + cells[nD].death + cells[nB].death
        sumant = sumant + cells[nD].epnumber + cells[nB].epnumber
        topant = maximum([antValues[nD], antValues[nB], topant])

        push!(tvec, t)
        push!(totAnt, sumant)
        push!(maxAnt, topant)

    end
    
    println("Simulation done.")
    return tvec, cells, totAnt, maxAnt, muts
end

function process_mutations(cells, detLim)
    antVec = Float64[]
    mutVec = Int64[]
    escVec = Bool[]
    for i=1:length(cells)
        append!(antVec, cells[i].epnumber)
        append!(mutVec, cells[i].mutations)
        append!(escVec, cells[i].escaped)
    end

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))

    return antVec, escVec, detMutDict

end

for i=1:20
    tvec, cells, totAnt, maxAnt, neoep_muts = birthdeath_neoep_constant(maxIter, popSize, p,pesc,neoep_dist, initialAnt, mu);
    outNDFsim = DataFrame(t=tvec, popAnt=totAnt, maxAnt=maxAnt)
    writetable("time_course_"*string(i)*".txt", outNDFsim)

    writedlm("neoep_mutations_"*string(i)*".txt", neoep_muts) #List of neoantigen mutations
    antVec, escVec, detMutDict = process_mutations(cells, detLim)
    writedlm("vaf_final_"*string(i)*".txt",detMutDict)
    popDF = DataFrame(antigenicity=antVec, escaped=escVec)
	writetable("final_population_"*string(i)*".txt",popDF)

end

