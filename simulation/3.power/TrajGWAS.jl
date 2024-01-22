# cd "/gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/"
# sbatch -J TrajGWAS --mem=2G -t 1-0:0 --array=1-1200 -o log/%A_%a.log --wrap='julia /gdata01/user/xuhe/family_relatedness/simulation-2023-09-19/code/power_TrajGWAS-2023-11-27.jl $SLURM_ARRAY_TASK_ID'

using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, NLopt, WiSER
using LinearAlgebra
using Serialization
using KNITRO

BLAS.set_num_threads(1)

n_cpu = parse(Int, ARGS[1])
table = DataFrame(rep = repeat(1:200, inner = 1, outer = 6),
 type = repeat(["A", "B", "C", "D", "E", "F"], inner = 200),
 scr = repeat(3:3, outer = 1200))
nrep = table[n_cpu, "rep"]; type = table[n_cpu, "type"]; scr = table[n_cpu, "scr"]
println("n.cpu = $n_cpu ; n.rep = $nrep ; type is $type ; scr is $scr.")

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

nulldf = CSV.read("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario$scr/pheno_new/pheno$type$scr-$nrep.csv", DataFrame)

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
f = open("/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario$scr/residuals_new/Resid$type-$nrep.txt", "w")
writedlm(f, ["R_beta" "R_tau" "SubjID"])

for j in 1:length(nm.data)
        R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
        R_tau = - sum(nm.data[j].diagDVRV)
        writedlm(f, [R_beta R_tau rownames[j]])
end
close(f)

filter!(:UNRELATED => UNRELATED -> UNRELATED == true, nulldf)

@time trajgwas(@formula(pheno ~ 1 + xone + xtwo + xthree),
                    @formula(pheno ~ 1 + zone),
                    @formula(pheno ~ 1 + xone + xtwo + xthree),
                    :SubjID,
                    nulldf,
                    "/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/GenoMat$type";
                    pvalfile="/gdata01/user/xuhe/SPA-GRM/simulation-2023-09-19/scenario$scr/TrajGWASpval_new/powerpval$type-$nrep.txt",
                    solver=solver,
                    solver_config = solver_config,
                    runs=10,
                    usespa=true)

println("completed")
