
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/"
# sbatch -J PRS --mem=15G -t 3-0:0 -o log/test_for_PRS_longitudinal_ferritin.log --wrap='julia /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/4.PRS_new/select_markers_longitudinal_new_step4_ferritin.jl'
using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using BGEN
using Serialization
using KNITRO

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

pheno = "ferritin"
println(pheno)

pheno_data = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/$pheno.txt", DataFrame)

pheno_data = filter(x -> x.BMI !== missing, pheno_data)

# standardize
standardizes(x) = (x .- mean(skipmissing(x))) ./ std(skipmissing(x))

pheno_data[!, :std_sf] = standardizes(pheno_data[!, :sf])
pheno_data[!, :std_age] = standardizes(pheno_data[!, :age])
pheno_data[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(pheno_data))
pheno_data[!, :std_age_sex] = pheno_data[!, :std_age] .* pheno_data[!, :sex_genetic]
pheno_data[!, :std_bmi] = standardizes(pheno_data[!, :BMI])

pheno_data[!, :age_indicator] = map(x -> round(x.age), eachrow(pheno_data))
pheno_data[!, :std_age_indicator] = standardizes(pheno_data[!, :age_indicator])

for chr in 1:22
    println(chr)

    LOCO_PRS = CSV.read(join(["/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/PRS/", pheno, chr, ".LOCO.txt"]), DataFrame)

    pheno_data_PRS = innerjoin(pheno_data, LOCO_PRS, on = :eid => :FID)

    # use WiSER to fit null model, note that we call WiSER by trajgwas without given genofile.
    @time nm = trajgwas(@formula(std_sf ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi + PRS +
            PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
        # @formula(std_sf ~ 1 + std_age_indicator),
        @formula(std_sf ~ 1),
        @formula(std_sf ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi),
        :eid,
        pheno_data_PRS,
        nothing;
        nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/$pheno.$chr.nullmodel.txt",
        solver=solver,
        solver_config = solver_config,
        runs=30,
        init = x -> WiSER.init_ls!(x; gniters=0))

    # extract the residuals for beta from fitted null model, note that the residuals must match id!!!
    SubjID = unique(nm.ids)
    f1 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/$pheno.$chr.betaresid.txt", "w")
    writedlm(f1, ["SubjID" "Resid"])

    for j in 1:length(nm.data)
            R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
            writedlm(f1, [SubjID[j] R_beta])
    end

    close(f1)

    # extract the residuals for tau from fitted null model, note that the residuals must match id!!!
    f2 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/ResidMat/$pheno.$chr.tauresid.txt", "w")
    writedlm(f2, ["SubjID" "Resid"])

    for j in 1:length(nm.data)
            R_tau = - sum(nm.data[j].diagDVRV)
            writedlm(f2, [SubjID[j]  R_tau])
    end

    close(f2)
end
