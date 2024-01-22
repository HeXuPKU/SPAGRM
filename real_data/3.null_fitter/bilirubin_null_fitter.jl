
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
solver = Ipopt.Optimizer(); solver_config = Dict("print_level"=>0, "watchdog_shortened_iter_trigger"=>3, "max_iter"=>100)
# solver = KNITRO.Optimizer(); solver_config = Dict("outlev"=>3) # (Knitro is commercial software)
# solver = NLopt.Optimizer(); solver_config = Dict("algorithm"=>:LD_MMA, "maxeval"=>4000)
# solver = NLopt.Optimizer(); solver_config = Dict("algorithm"=>:LD_LBFGS, "maxeval"=>4000)

bilirubin_data = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bilirubin.txt", DataFrame)

bilirubin_data = filter(x -> x.BMI !== missing, bilirubin_data)

# standardize
standardizes(x) = (x .- mean(skipmissing(x))) ./ std(skipmissing(x))

bilirubin_data[!, :std_bili] = standardizes(bilirubin_data[!, :bili])
bilirubin_data[!, :std_age] = standardizes(bilirubin_data[!, :age])
bilirubin_data[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(bilirubin_data))
bilirubin_data[!, :std_age_sex] = bilirubin_data[!, :std_age] .* bilirubin_data[!, :sex_genetic]
bilirubin_data[!, :std_bmi] = standardizes(bilirubin_data[!, :BMI])

bilirubin_data[!, :age_indicator] = map(x -> round(x.age), eachrow(bilirubin_data))
bilirubin_data[!, :std_age_indicator] = standardizes(bilirubin_data[!, :age_indicator])

# use WiSER to fit null model, note that we call WiSER by trajgwas without given genofile.
@time nm = trajgwas(@formula(std_bili ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_bili ~ 1 + std_age_indicator),
    @formula(std_bili ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi),
    :eid,
    bilirubin_data,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_nm.txt",
    solver=solver,
    solver_config = solver_config,
    runs=10,
    init = x -> WiSER.init_ls!(x; gniters=0))

# extract the residuals for beta from fitted null model, note that the residuals must match id!!!
SubjID = unique(nm.ids)
f1 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_beta_resid.txt", "w")
writedlm(f1, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
        writedlm(f1, [SubjID[j] R_beta])
end

close(f1)

# extract the residuals for tau from fitted null model, note that the residuals must match id!!!
f2 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_tau_resid.txt", "w")
writedlm(f2, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_tau = - sum(nm.data[j].diagDVRV)
        writedlm(f2, [SubjID[j]  R_tau])
end

close(f2)

# extract unrelated subjects to run trajgwas for latter comparison.
bilirubin_data_un = filter(x -> x.UNRELATED_WB == true, bilirubin_data)

bilirubin_data_un[!, :std_bili] = standardizes(bilirubin_data_un[!, :bili])
bilirubin_data_un[!, :std_age] = standardizes(bilirubin_data_un[!, :age])
bilirubin_data_un[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(bilirubin_data_un))
bilirubin_data_un[!, :std_age_sex] = bilirubin_data_un[!, :std_age] .* bilirubin_data_un[!, :sex_genetic]
bilirubin_data_un[!, :std_bmi] = standardizes(bilirubin_data_un[!, :BMI])

bilirubin_data_un[!, :age_indicator] = map(x -> round(x.age), eachrow(bilirubin_data_un))
bilirubin_data_un[!, :std_age_indicator] = standardizes(bilirubin_data_un[!, :age_indicator])

@time nm_un = trajgwas(@formula(std_bili ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi + 
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_bili ~ 1 + std_age_indicator),
    @formula(std_bili ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi),
    :eid,
    bilirubin_data_un,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_nm_un.txt",
    solver=solver,
    solver_config = solver_config,
    runs=10,
    init = x -> WiSER.init_ls!(x; gniters=0))

open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_nm_un.jls", "w") do io
        Serialization.serialize(io, nm_un)
end

# using Serialization
# nm_un = open(deserialize, "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/bilirubin_nm_un.jls")