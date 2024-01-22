
# cd "/gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/"
# sbatch -J TrajGWAS --mem=5G -t 1-0:0 --array=1-22 -o log/%A_%a.log --wrap='julia /gdata01/user/xuhe/family_relatedness/real_data-2022-12-29/4.trajgwas/FVC_TrajGWAS.jl $SLURM_ARRAY_TASK_ID'
n_cpu = parse(Int, ARGS[1])
print(n_cpu)

using DataFrames, CSV, Dates, DelimitedFiles
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using BGEN
using Serialization

### GWAS for whole genome.
nm = open(deserialize, "/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/null_fitter/FVC_nm_un.jls")
genetic_iids_subsample = nm.ids

bgenfile = join(["/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c", n_cpu, "_b0_v3.bgen"])
samplefile = "/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb22828_c1_b0_v3_s487203.sample"
mfifile = join(["/gdata02/master_data1/UK_Biobank/ukb22828_imp/ukb_mfi_chr", n_cpu, "_v3.txt"])

pvalfile = join(["/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/temp_pvalues/FVC/chr", n_cpu, ".trajgwas.txt"])

ukb_data = Bgen(bgenfile; sample_path = samplefile)
genetic_iids = map(x -> parse(Int, split(x, " ")[1]), samples(ukb_data))

order_dict = Dict{Int, Int}()
for (i, iid) in enumerate(genetic_iids)
    order_dict[iid] = i
end

sample_indicator = falses(length(genetic_iids))
for v in genetic_iids_subsample
    sample_indicator[order_dict[v]] = true # extract only the samples being used for the analysis
end

mfi = CSV.read(mfifile, DataFrame, header = false)
mfi.Column8 = map(x -> x == "NA" ? NaN : parse(Float64, x), mfi.Column8)
snpmask = (mfi.Column6 .> 0.0005) .& (mfi.Column8 .> 0.6)

# rearrange data in nm so that it matches bgen data
nullinds = indexin(genetic_iids[sample_indicator], nm.ids)
nm.obswts .= isempty(nm.obswts) ? nm.obswts : nm.obswts[nullinds]
nm.ids .= nm.ids[nullinds]
nm.nis .= nm.nis[nullinds]
nm.data .= nm.data[nullinds]
@assert genetic_iids[sample_indicator] == nm.ids "there is some issue -- sampleids not matching"

@time trajgwas(nm, 
    bgenfile, 
    count(sample_indicator);
    samplepath = samplefile,
    pvalfile = pvalfile,
    snpinds = snpmask,
    bgenrowinds = sample_indicator,
    usespa = true)
