# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/"
# sbatch -J TrajGWAS --mem=5G -t 1-0:0 --array=1-30000 -o log/%A_%a.log --wrap='julia /gdata01/user/xuhe/family_relatedness/simulation-2023-03-22/code/type1error_TrajGWAS-2023-07-01.jl $SLURM_ARRAY_TASK_ID'

using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using Serialization
using KNITRO

n_cpu = parse(Int, ARGS[1])
table = DataFrame(rep = repeat(1:10000, outer = 3),
 type = repeat(["A", "B", "C"], inner = 10000),
 genofile = repeat(["nSub_25000_nFam_6250", "nSub_25000_nFam_2500", "nSub_50000"], inner = 10000))
nrep = table[n_cpu, "rep"]; type = table[n_cpu, "type"]; genofile = table[n_cpu, "genofile"]
println("n.cpu = $n_cpu ; n.rep = $nrep ; type is $type; genofile is $genofile.")

BLAS.set_num_threads(1)

# for TrajGWAS old version, solvers settings:
# solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes", warm_start_init_point="yes", max_iter=100)
# solver = Ipopt.IpoptSolver(print_level=0, watchdog_shortened_iter_trigger=3, max_iter=100)
# solver = KNITRO.KnitroSolver(outlev=3) # (Knitro is commercial software)
# solver = NLopt.NLoptSolver(algorithm=:LD_MMA, maxeval=4000)
# solver = NLopt.NLoptSolver(algorithm=:LD_LBFGS, maxeval=4000)

# for TrajGWAS new version, solvers settings:
# solver = Ipopt.Optimizer(); solver_config = Dict("print_level"=>0, "mehrotra_algorithm"=>"yes", "warm_start_init_point"=>"yes", "max_iter"=>100)
# solver = Ipopt.Optimizer(); solver_config = Dict("print_level"=>0, "watchdog_shortened_iter_trigger"=>3, "max_iter"=>100)
# solver = KNITRO.Optimizer(); solver_config = Dict("outlev"=>3) # (Knitro is commercial software)
solver = NLopt.Optimizer(); solver_config = Dict("algorithm"=>:LD_MMA, "maxeval"=>4000)
# solver = NLopt.Optimizer(); solver_config = Dict("algorithm"=>:LD_LBFGS, "maxeval"=>4000)

nulldf = CSV.read("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/pheno/pheno$type-$nrep.csv", DataFrame)

@time nm = trajgwas(@formula(pheno ~ 1 + xone + xtwo + xthree),
                   @formula(pheno ~ 1 + zone),
                   @formula(pheno ~ 1 + xone + xtwo + xthree),
                   :SubjID,
                   nulldf,
                   nothing;
                   solver=solver,
                   solver_config = solver_config,
                   runs=10)

println(nm)

rownames = unique(nm.ids)
f = open("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/residuals/Resid$type-$nrep.txt", "w")
writedlm(f, ["R_beta" "R_tau" "SubjID"])

for j in 1:length(nm.data)
        R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
        R_tau = - sum(nm.data[j].diagDVRV)
        writedlm(f, [R_beta R_tau rownames[j]])
end
close(f)

trajgwas(nm,
        join(["/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/", genofile, "common/", genofile]);
        pvalfile = "/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/TrajGWASpval/commonpval$type-$nrep.txt",
        usespa = true)

trajgwas(nm,
        join(["/gdata01/user/xuhe/SPA-GRM/SimulatedGenotypeUsingWES/", genofile, "rare/", genofile]);
        pvalfile = "/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/TrajGWASpval/rarepval$type-$nrep.txt",
        usespa = true)

# rm("/gdata01/user/xuhe/SPA-GRM/simulation-2023-03-22/scenario1/pheno/pheno$type-$nrep.csv")

println("completed")
