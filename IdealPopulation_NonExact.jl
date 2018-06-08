using Distributions, StatsBase, DataFrames, GLM

type cancercell
    mutations::Array{Int64,1}
    fitness::Int64
end

function newmutations(cancercell, mutID, p)

    numbermutations= 1
    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    
    neoep = rand()<p
    neoep_effective = neoep * (rand()<0.1)
    if neoep_effective
        cancercell.fitness = 0
    end

    return cancercell, mutID, neoep
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness))
end

function tumourgrow_birthdeath_neoep(b0, d0, b_, d_, Nmax, p)


    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b0+d_

    #initialize arrays and parameters
    mutID = 1
    cells = cancercell[]
    muts = Int64[]
    neoep_muts = Int64[]
    push!(cells,cancercell([],1))
    for i=1:10
    cells[1],mutID,neoep = newmutations(cells[1],mutID, p)
    #cells[1].fitness = 1 #manually overwrite
    push!(muts,mutID)
    end
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        if cells[randcell].fitness==0 #set death rate according to whether cell is antigenic or not
            d = d_
        else
            d = d0
        end

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #add new mutations to both new cells
            cells[randcell],mutID,neoep = newmutations(cells[randcell],mutID, p)
            push!(muts,mutID)
            if neoep
                push!(neoep_muts, mutID)
            end
            cells[end],mutID,neoep = newmutations(cells[end],mutID, p)
            push!(muts,mutID)
            if neoep
                push!(neoep_muts, mutID)
            end
            
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        if  (b0+d)<= r
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if b0 <= r < (b0+d)

            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies reinitialize simulation
        if (N == 0)
                mutID = 1
        cells = cancercell[]
        muts = Int64[]
    neoep_muts = Int64[]
    push!(cells,cancercell([],1))
            for i=1:10
    cells[1],mutID,neoep = newmutations(cells[1],mutID, p)
    #cells[1].fitness = 1 #manually overwrite
    push!(muts,mutID)
            end
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
        end

    end
    
    return Nvec, tvec, muts, neoep_muts, cells
end

function process_mutations(cells, neoep_muts, detLim, popSize, mu)
	fitnessVec = []
	mutVec = []
    for i=1:length(cells)
        push!(fitnessVec, cells[i].fitness)
        append!(mutVec, cells[i].mutations)
    end
	println("Non-immunogenic cells: ", sum(fitnessVec))

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))
    detEpMutDict = filter((k, v) -> k in neoep_muts, detMutDict)
    VAF = collect(values(detMutDict))/(2*popSize)
    VAFtotal = Array{Float64}(0)
	for mut in VAF
  	  n = rand(Poisson(mu))
  	  append!(VAFtotal, repeat([mut], inner=n))
	end
    append!(VAFtotal, VAF)
    VAFep = collect(values(detEpMutDict))/(2*popSize)
	println("Number of neo-epitope mutations: ",length(detEpMutDict))
	return length(detEpMutDict)/length(detMutDict), VAFtotal, VAFep

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
outREpvec = []
outMutRvec = []
outvafDF = DataFrame(invf = (1./steps - 1/fmax))
outvafEpDF = DataFrame(invf = (1./steps - 1/fmax))
outfimDF = DataFrame(f=steps)
outfimEpDF = DataFrame(f=steps)
for i=1:200
    Nvec, tvec, muts, neoep_muts, cells = tumourgrow_birthdeath_neoep(1, d0, 1, d_, popSize, p);
    mutR, VAFtotal, VAFep = process_mutations(cells, neoep_muts, detLim, popSize, mu)
    cumcount, cumfm, rsq = analyse_vaf(VAFtotal, steps, fmin, fmax)
    cumcountEp, cumfmEp, rsqEp = analyse_vaf(VAFep, steps, fmin, fmax)
    push!(outMutRvec, mutR)
    push!(outRvec, rsq)
    push!(outREpvec, rsqEp)
    outvafDF[Symbol("sim_$(i)")] = cumcount
    outfimDF[Symbol("sim_$(i)")] = cumfm
    outvafEpDF[Symbol("sim_$(i)")] = cumcountEp
    outfimEpDF[Symbol("sim_$(i)")] = cumfmEp
end

writetable("Vafdf.csv",outvafDF)
writetable("Fimdf.csv",outfimDF)
writetable("Vafdf_ep.csv",outvafEpDF)
writetable("Fimdf_ep.csv",outfimEpDF)
writedlm("Rsqs.txt",outRvec)
writedlm("Rsqs_ep.txt", outREpvec)
writedlm("MutRatios.txt", outMutRvec)
