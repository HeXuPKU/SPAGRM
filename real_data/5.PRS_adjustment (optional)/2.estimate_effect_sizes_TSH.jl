
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/"
# sbatch -J PRS --mem=15G -t 3-0:0 -o log/test_for_PRS_longitudinal_TSH.log --wrap='julia /gdata01/user/xuhe/family_relatedness/real_data-2023-08-29/4.PRS_new/select_markers_longitudinal_new_step2_TSH.jl'
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

pheno = "TSH"
println(pheno)

pheno_data = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/$pheno.txt", DataFrame)

pheno_data = filter(x -> x.BMI !== missing, pheno_data)
pheno_data = filter(x -> x.UNRELATED_WB == true, pheno_data)
genetic_iids_subsample = unique(pheno_data.eid)

# standardize
standardizes(x) = (x .- mean(skipmissing(x))) ./ std(skipmissing(x))

pheno_data[!, :std_tsh] = standardizes(pheno_data[!, :tsh])
pheno_data[!, :std_age] = standardizes(pheno_data[!, :age])
pheno_data[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(pheno_data))
pheno_data[!, :std_age_sex] = pheno_data[!, :std_age] .* pheno_data[!, :sex_genetic]
pheno_data[!, :std_bmi] = standardizes(pheno_data[!, :BMI])

pheno_data[!, :age_indicator] = map(x -> round(x.age), eachrow(pheno_data))
pheno_data[!, :std_age_indicator] = standardizes(pheno_data[!, :age_indicator])

for chr in 1:22
    println(chr)

    bgenfile = join(["/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", chr, "_b0_v3"])
    samplefile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"

    if(!isfile("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/$pheno.$chr.snpmask.txt"))
        continue
    end

    snpmask = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/$pheno.$chr.snpmask.txt", DataFrame, header = false)
    snpmask = snpmask.Column1

    ukb_data = Bgen(bgenfile * ".bgen"; sample_path = samplefile)
    genetic_iids = map(x -> parse(Int, split(x, " ")[1]), samples(ukb_data))

    order_dict = Dict{Int, Int}()
    for (i, iid) in enumerate(genetic_iids)
        order_dict[iid] = i
    end

    sample_indicator = falses(length(genetic_iids))
    for v in genetic_iids_subsample
        sample_indicator[order_dict[v]] = true # extract only the samples being used for the analysis
    end

    @time trajgwas(@formula(std_tsh ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    # @formula(std_tsh ~ 1 + std_age_indicator),
    @formula(std_tsh ~ 1), # more efficient
    @formula(std_tsh ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi),
    :eid,
    pheno_data,
    bgenfile;
    samplepath = samplefile,
    geneticformat = "BGEN",
    pvalfile = "/gdata01/user/xuhe/SPA-GRM/real_data-2023-08-29/PRS_new/longitudinal/lead_SNPs/$pheno.$chr.effectsize.txt",
    snpinds = snpmask,
    geneticrowinds = sample_indicator,
    test = :wald,
    solver=solver,
    solver_config = solver_config,
    runs=10,
    init = x -> WiSER.init_ls!(x; gniters=0),
    disable_wsvar = true)
end
