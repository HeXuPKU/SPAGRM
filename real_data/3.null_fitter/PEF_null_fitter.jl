
using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
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

PEF_data = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/PEF.txt", DataFrame)

PEF_data = filter(x -> x.BMI !== missing, PEF_data)
PEF_data = filter(x -> x.genoArray !== missing, PEF_data)
PEF_data = filter(x -> x.smoke !== missing, PEF_data)

# standardize
standardizes(x) = (x .- mean(skipmissing(x))) ./ std(skipmissing(x))

PEF_data[!, :std_pef] = standardizes(PEF_data[!, :pef])
PEF_data[!, :std_age] = standardizes(PEF_data[!, :age])
PEF_data[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(PEF_data))
PEF_data[!, :std_age_sex] = PEF_data[!, :std_age] .* PEF_data[!, :sex_genetic]
PEF_data[!, :std_bmi] = standardizes(PEF_data[!, :BMI])

PEF_data[!, :age_indicator] = map(x -> round(x.age), eachrow(PEF_data))
PEF_data[!, :std_age_indicator] = standardizes(PEF_data[!, :age_indicator])

# use WiSER to fit null model, note that we call WiSER by trajgwas without given genofile.
@time nm = trajgwas(@formula(std_pef ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        genoArray + smoke + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    # @formula(std_pef ~ 1 + std_age_indicator),
    @formula(std_pef ~ 1),
    @formula(std_pef ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi + genoArray + smoke),
    :eid,
    PEF_data,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_nm.txt",
    solver=solver,
    solver_config = solver_config,
    runs=10)

# extract the residuals for beta from fitted null model, note that the residuals must match id!!!
SubjID = unique(nm.ids)
f1 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_beta_resid.txt", "w")
writedlm(f1, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
        writedlm(f1, [SubjID[j] R_beta])
end

close(f1)

# extract the residuals for tau from fitted null model, note that the residuals must match id!!!
f2 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_tau_resid.txt", "w")
writedlm(f2, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_tau = - sum(nm.data[j].diagDVRV)
        writedlm(f2, [SubjID[j]  R_tau])
end

close(f2)

# extract unrelated subjects to run trajgwas for latter comparison.
PEF_data_un = filter(x -> x.UNRELATED_WB == true, PEF_data)

PEF_data_un[!, :std_pef] = standardizes(PEF_data_un[!, :pef])
PEF_data_un[!, :std_age] = standardizes(PEF_data_un[!, :age])
PEF_data_un[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(PEF_data_un))
PEF_data_un[!, :std_age_sex] = PEF_data_un[!, :std_age] .* PEF_data_un[!, :sex_genetic]
PEF_data_un[!, :std_bmi] = standardizes(PEF_data_un[!, :BMI])

PEF_data_un[!, :age_indicator] = map(x -> round(x.age), eachrow(PEF_data_un))
PEF_data_un[!, :std_age_indicator] = standardizes(PEF_data_un[!, :age_indicator])

@time nm_un = trajgwas(@formula(std_pef ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        genoArray + smoke + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    # @formula(std_pef ~ 1 + std_age_indicator),
    @formula(std_pef ~ 1),
    @formula(std_pef ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi + genoArray + smoke),
    :eid,
    PEF_data_un,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_nm_un.txt",
    solver=solver,
    solver_config = solver_config,
    runs=10)

open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_nm_un.jls", "w") do io
        Serialization.serialize(io, nm_un)
end

# using Serialization
# nm_un = open(deserialize, "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/PEF_nm_un.jls")