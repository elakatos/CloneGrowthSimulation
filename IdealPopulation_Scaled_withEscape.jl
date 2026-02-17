#=
    Julia code implementing simulations of a stochastic branching process model of tumour growth and mutation accumulations, under negative selection from the immune system due to randomly arising antigenic mutations. This code corresponds to an ideal population/measurement, with no measurement noise (but intrinsic stochasticity) and user-defined detection limit for mutations.

    For details of the model see Lakatos et al, Nat Genet, 2020: 10.1038/s41588-020-0687-1 and the Readme at https://github.com/elakatos/CloneGrowthSimulation

    This file implements the base simulations with probabilistic active immune escape, acquired by chance throughout tumour growth.
    Parameters: d0, popSize, detLim, p, initial_mut, mu, neoep_dist, pesc
    Outputs:
        - mutations_<i>.txt : data frame of mutations with their antigenicity value, escape status (true/false) and number of cells the mutation was observed. Mutations below detection limit are outputted
        - population_<i>.txt : cell growth curve with time, total population and non-immunogenic population in three comma-separated columns

    @Author: Eszter Lakatos (eszter.lakatos@chalmers.se)
=#


using ArgParse
using Distributions
using StatsBase
using Random
using CSV
using DelimitedFiles
using DataFrames

# ==============================
# Parameter container and initialisation
# ==============================
struct SimParams{F}
    b0::Float64
    d0::Float64
    Nmax::Int
    p::Float64
    pesc::Float64
    mu::Float64
    s::Float64
    fitness_function::F
    immThresh::Float64
    tmax::Float64
    detLim::Int
    imut::Int
    clonalEscape::Bool
    clonalImm::Float64
    nSim::Int
    folder::String
end

function parse_parameters()

    s_params = ArgParseSettings()

    @add_arg_table s_params begin
        "--s"
            arg_type = Float64
            required = true
            help = "Selection against immunogenic cells (negative value)"
       "--Nmax"
            arg_type = Int
            default = 1e5
            help = "Maximum cancer population size"
        "--mu"
            arg_type = Float64
            default = 1
            help = "Mutation rate"
        "--imut"
            arg_type = Int
            default = 10
            help = "Number of initial mutations in founder cell"
        "--p"
            arg_type = Float64
            default = 0.075
            help = "Probability of immunogenic mutation"
        "--pesc"
            arg_type = Float64
            default = 1e-6
            help = "Probability of escape mutation"
        "--fitness"
            arg_type = String
            default = "linear"
            help = "Name of fitness function to use, available options are: linear / threshold / log"
        "--b0"
            arg_type = Float64
            default = 1.0
            help = "Birth rate"
        "--d0"
            arg_type = Float64
            default = 0.1
            help = "Baseline (non-immunogenic) death rate"
        "--immThresh"
            arg_type = Float64
            default = 0.5
            help = "Threshold of a cell being immunogenic (for visualisation)"
        "--tmax"
            arg_type = Float64
            default = 300
            help = "Maximum simulation time"
        "--detLim"
            arg_type = Int
            default = 5
            help = "Lower limit of mutation detection, in number of cells carrying mutation"
        "--cEscape"
            arg_type = Bool
            default = false
            help = "Boolean to deterministically set if the founder cell is escaped"
        "--cImm"
            arg_type = Float64
            default = -1.0
            help = "Deterministically set value of the founder cell's immunogenicity; negative values keep the assignment random"
        "--out"
            arg_type = String
            default = "."
            help = "Output folder"
        "--nSim"
            arg_type = Int
            default = 100
            help = "Number of simulated tumours to create"
    end

    parsed = parse_args(s_params)

    return SimParams(
        parsed["b0"],
        parsed["d0"],
        parsed["Nmax"],
        parsed["p"],
        parsed["pesc"],
        parsed["mu"],
        parsed["s"],
        get_fitness_function(parsed["fitness"]),
        parsed["immThresh"],
        parsed["tmax"],
        parsed["detLim"],
        parsed["imut"],
        parsed["cEscape"],
        parsed["cImm"],
        parsed["nSim"],
        parsed["out"]
    )
end

# ==============================
# Fitness function
# ==============================

linear_fitness(n, s) = 1 + s*n
threshold_fitness(n, s) = n > 1.0 ? 0.0 : 1.0
log_fitness(n, s) = 1 + s*log(2*n+1)

function get_fitness_function(name)
    if name == "linear"
        return linear_fitness
    elseif name == "log"
        return log_fitness
    elseif name == "threshold"
        return threshold_fitness
    else
        error("Unknown fitness function: $name")
    end
end

# ==============================
# Mutation storage
# ==============================
mutable struct MutationStore
    immune::Vector{Float64}
    escape::Vector{Bool}
    count::Vector{Int}
end

# ==============================
# Cancer cell
# ==============================
mutable struct CancerCell
    epnumber::Float64
    escaped::Bool
    nonimmunogenic::Bool
    mutations::Vector{Int}
end

# ==============================
# Add mutations (in-place)
# ==============================
function add_mutations!(
    cell::CancerCell,
    store::MutationStore,
    mutID::Int,
    params::SimParams,
    neoep_dist,
    n_new::Int
)
    @inbounds for _ in 1:n_new
        push!(cell.mutations, mutID)
        push!(store.count, 1)

        # Neoepitope mutation
        if rand() < params.p
            neoep_val = max(0.0, rand(neoep_dist))
            push!(store.immune, neoep_val)
            cell.epnumber += neoep_val
            cell.nonimmunogenic = cell.epnumber < params.immThresh
        else
            push!(store.immune, 0.0)
        end

        # Escape mutation
        if rand() < params.pesc
            cell.escaped = true
            push!(store.escape, true)
        else
            push!(store.escape, false)
        end
        mutID += 1
    end

    return mutID
end

# ==============================
# Initialize population
# ==============================
function start_population(params::SimParams, neoep_dist)

    cells = CancerCell[]
    sizehint!(cells, params.Nmax)

    mut_store = MutationStore(Float64[],Bool[],Int[])
    mutID = 1

    if params.clonalImm >= 0.0 # if there is a pre-defined immunogenicity value for the first cell
        if params.clonalEscape
            cell = CancerCell(params.clonalImm, true, params.clonalImm < params.immThresh, Int[])
        else
            cell = CancerCell(params.clonalImm, false, params.clonalImm < params.immThresh, Int[])
        end
    else # otherwise, add mutations and immunogenicity (and potentially escape) is decided stochastically
        if params.clonalEscape
            cell = CancerCell(0.0, true, true, Int[])
        else
            cell = CancerCell(0.0, false, true, Int[])
        end
        mutID = add_mutations!(
            cell, mut_store, mutID, params, neoep_dist, params.imut
        )

    end

    push!(cells, cell)
    nonimm = 0
    nonimm = cell.nonimmunogenic

    return cells, mutID, mut_store, nonimm

end

# ==============================
# Main simulation
# ==============================
function birthdeath_neoep(
    params::SimParams,
    neoep_dist;
    initial_mut=10,
)

    cells, mutID, mut_store, nonimm = 
        start_population(params, neoep_dist)

    pois = Poisson(params.mu)

    N = 1
    t = 0.0
    dmax = params.d0

    Nvec = Int[]
    tvec = Float64[]
    nonimmvec = Int[]

    # allocation hint: expected number of steps increases with Nmax and is typically 5-10x that
    sizehint!(Nvec, 10*params.Nmax)
    sizehint!(tvec, 10*params.Nmax)
    sizehint!(nonimmvec, 10*params.Nmax)

    push!(Nvec, N)
    push!(tvec, t)
    push!(nonimmvec, nonimm)

    while (N < params.Nmax) && (t < params.tmax)
        # choose a random cell that will undergo some event
        randcell = rand(1:N)
        cell = cells[randcell]
        Nt = N

        # death rate, defined based on escape status and fitness (note that fitness can even be <0)
        cell_fitness = params.fitness_function(cell.epnumber, params.s)
        d = cell.escaped ? params.d0 :
            max(0.0, (params.d0 - params.b0) * cell_fitness + params.b0)

        dmax = max(dmax, d)
        Rmax = params.b0 + dmax

        r = rand() * Rmax

        # ================= Birth =================
        if r < params.b0

            N += 1

            # copy
            newcell = CancerCell(
                cell.epnumber,
                cell.escaped,
                cell.nonimmunogenic,
                copy(cell.mutations),
            )
            push!(cells, newcell)
            for m in cell.mutations # increase count of how many cells carry each mutation
                mut_store.count[m] += 1
            end

            nonimm -= cell.nonimmunogenic

            mutID = add_mutations!(
                cell, mut_store, mutID, params, neoep_dist, rand(pois)
            )
            mutID = add_mutations!(
                newcell, mut_store, mutID, params, neoep_dist, rand(pois)
            )

            nonimm += cell.nonimmunogenic + newcell.nonimmunogenic
        end

        # ================= Death =================
        if (params.b0 <= r) && (r < params.b0 + d)

            N -= 1
            nonimm -= cell.nonimmunogenic

            for m in cell.mutations #decrease count of how many cells carry each mutation
                mut_store.count[m] -= 1
            end

            # removal with pop
            cells[randcell] = cells[end]
            pop!(cells)
        end

        # ================= Time step =================
        t += randexp() / (Rmax * Nt)

        push!(Nvec, N)
        push!(tvec, t)
        push!(nonimmvec, nonimm)

        # restart if extinct
        if N == 0
            cells, mutID, mut_store, nonimm =
                start_population(params, neoep_dist)
            N = 1
            dmax = params.d0     # when population restarts, we reset dmax as well
            push!(Nvec, N)
            push!(tvec, t)
            push!(nonimmvec, nonimm)
        end
    end

    return Nvec, tvec, mutID, mut_store, cells, nonimmvec

end


params = parse_parameters()
println("Running simulation with:")
println(params)

neoep_dist = Exponential(0.2)

print("Simulation number: ")
for i = 1:params.nSim
    Nvec, tvec, mutID, mut_store, cells, immune = birthdeath_neoep(params, neoep_dist)
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    CSV.write(params.folder*"/population_"*string(i)*".csv", outNDFsim)
    detected_muts = filter(i -> mut_store.count[i] > params.detLim, eachindex(mut_store.count)) # filter for mutations above detection limit
    outMutDFsim = DataFrame(imm=mut_store.immune[detected_muts], esc=mut_store.escape[detected_muts], count=mut_store.count[detected_muts])
    CSV.write(params.folder*"/mutations_"*string(i)*".csv", outMutDFsim)
    print(string(i)*"..")
    epnums = [c.epnumber for c in cells]
    writedlm(params.folder*"/cell_immunogenicities_"*string(i)*".txt", epnums)
end

println("\nFinished simulation and saved simulation results.")
