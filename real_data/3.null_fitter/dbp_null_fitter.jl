
using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using Serialization
# using KNITRO

BLAS.set_num_threads(1)
# solver = Ipopt.IpoptSolver(print_level=0, mehrotra_algorithm = "yes", warm_start_init_point="yes", max_iter=100)
# solver = Ipopt.Optimizer()
# solver = KNITRO.KnitroSolver(outlev=3) # (Knitro is commercial software)
solver = NLopt.NLoptSolver(algorithm=:LD_MMA, maxeval=4000)
# solver = NLopt.NLoptSolver(algorithm=:LD_LBFGS, maxeval=4000)
# solver = Ipopt.IpoptSolver(print_level=0, watchdog_shortened_iter_trigger=3, max_iter=100)

bp_data = CSV.read("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/bp.txt", DataFrame)

bp_data = filter(x -> x.BMI !== missing, bp_data)
bp_data = filter(x -> x.self_bpdrugs !== -1, bp_data)

# standardize
standardizes(x) = (x .- mean(skipmissing(x))) ./ std(skipmissing(x))

bp_data[!, :std_dbp] = standardizes(bp_data[!, :DBP_shifted])
bp_data[!, :std_age] = standardizes(bp_data[!, :age])
bp_data[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(bp_data))
bp_data[!, :std_age_sex] = bp_data[!, :std_age] .* bp_data[!, :sex_genetic]
bp_data[!, :std_bmi] = standardizes(bp_data[!, :BMI])

bp_data[!, :age_indicator] = map(x -> round(x.age), eachrow(bp_data))
bp_data[!, :std_age_indicator] = standardizes(bp_data[!, :age_indicator])

@time nm = trajgwas(@formula(std_dbp ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_dbp ~ 1 + std_age_indicator),
    @formula(std_dbp ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    :eid,
    bp_data,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_nm.txt",
    solver=solver,
    runs=10)

SubjID = unique(nm.ids)
# extract the residuals for beta from fitted null model, note that the residuals must match id!!!
f1 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_beta_resid.txt", "w")
writedlm(f1, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_beta = sum(nm.data[j].Dinv_r - transpose(nm.data[j].rt_UUt))
        writedlm(f1, [SubjID[j] R_beta])
end

close(f1)

# extract the residuals for tau from fitted null model, note that the residuals must match id!!!
f2 = open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_tau_resid.txt", "w")
writedlm(f2, ["SubjID" "Resid"])

for j in 1:length(nm.data)
        R_tau = - sum(nm.data[j].diagDVRV)
        writedlm(f2, [SubjID[j]  R_tau])
end

close(f2)


# extract unrelated subjects to run trajgwas for latter comparison.
bp_data_un = filter(x -> x.UNRELATED_WB == true, bp_data)

bp_data_un[!, :std_dbp] = standardizes(bp_data_un[!, :DBP_shifted])
bp_data_un[!, :std_age] = standardizes(bp_data_un[!, :age])
bp_data_un[!, :std_age_sq] = map(x -> x.std_age ^ 2, eachrow(bp_data_un))
bp_data_un[!, :std_age_sex] = bp_data_un[!, :std_age] .* bp_data_un[!, :sex_genetic]
bp_data_un[!, :std_bmi] = standardizes(bp_data_un[!, :BMI])

bp_data_un[!, :age_indicator] = map(x -> round(x.age), eachrow(bp_data_un))
bp_data_un[!, :std_age_indicator] = standardizes(bp_data_un[!, :age_indicator])

@time nm_un = trajgwas(@formula(std_dbp ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    @formula(std_dbp ~ 1 + std_age_indicator),
    @formula(std_dbp ~ 1 + sex_genetic + std_age + std_age_sq + std_age_sex + std_bmi +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10),
    :eid,
    bp_data_un,
    nothing;
    nullfile="/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_nm_un.txt",
    solver=solver,
    runs=10)

open("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_nm_un.jls", "w") do io
        Serialization.serialize(io, nm_un)
end

# using Serialization
# nm_un = open(deserialize, "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/dbp_nm_un.jls")
