
###### phecodes used in SPAGRM. 

# excluded_BP_codes <- '246[CDNPS]|^XaJ2[GH]|^XaIw[jk]' # containing standing and lying bp.
BP_codes <- "^246[.abcdefgABEFGJQRTVWXY12345679]|^XaF4[abFKLOSZ]|^XaJ2[EF]|^XaKF[xw]|^XaKj[FG]|^662L|^R1y2|^G20"

height_weight_BMI_codes <- '^XaCDR|^XaJJH|^XaJqk|^XaZcl|^22K|^229|^22A|^162[23]|^X76CG|^XE1h4|^XM01G|^Xa7wI' # no change compared to TrajGWAS.

waistcircu_codes <- '^22N0|^Xa041'

pulse_rate_codes <- '^242|^243|^X773s|^XaIBo|^XM02J'
CHDscore_codes <- '^38DP|^XaFsZ'

HDL_codes <- '^44d[23A]|^X772M|^44P[5BC]|^XaEVr' # omit 44R3 : HDL-protein
LDL_codes <- '^44d[45B]|^X772N|^44P[6DEI]|^XaEVs|^XaIp4' # omit 44R4 : HDL-protein
NHDL_codes <- '^44PL|^XabE1|^XaN3z'
totchol_codes <- '^44OE|^44P[.12349HJKZ]|^XE2eD|^X772L|^XSK14|^XaFs9|^XaIRd|^XaJe9|^XaLux'
triglyc_codes <- '^44e|^44Q|^X772O|^XE2q9' # no change compared to TrajGWAS.
# TC_HDL_codes <- '^44l[2FGH]|^44PF|^XaEU[qr]' # omit 44PG. and XaEil

fastgluc_codes <- '^44[fg]1|^44T[2K]|^XE2mq'
randgluc_codes <- '^44[fg][.0]|^44T[1AJ]|^44U\\.|^8A17|^X772z|^XE2mp|^XM0ly'
a1c_codes <- '^XaPbt|^XaERp|^X772q|^42W[.12345Z]\\.|^44TB\\.' # no change compared to TrajGWAS.

spo2_codes <- 'X770D'
FEV1_codes <- '^33972|^339[abef]|^339O\\.|^X77Qu|^XaIx[QRUV]'
FVC_codes <- '^3396[.013]|^339[hs]|^XaJ3K|^XaPpI'
FEV1_FVC_codes <- '^X77Ra|^XaEFy'
# excluded_PEF_codes <- '^745C0|^339[dB]|^XaEGA|^XaIxT|^XE2ws' # omit PEF after bronchodilation or steroids.
PEF_codes <- '^339[345Acgno]|^66YX|^X77RW|^XaEFM|^XaEHe|^XaIx[DPS]|^XaJEg|^XaK6I|^XE2wr|^XE2xQ|^XS7q5'

RBC_codes <- '^426[.123457Z]'
Platelet_codes <- '^42P[.1234Z]|^XE24o'
WBC_codes <- '^42H[.123578Z]|^XaId[YZ]'

BASO_codes <- '^42L.'
EOS_codes <- '^42K.'
NEUT_codes <- '^42J[.15Z]'
LYMPH_codes <- '^42M[.14Z]'
MONO_codes <- '^42N[.125Z]'
perc_BASO_codes <- '^42b3|^XE2mS'
perc_EOS_codes <- '^42b9|^XaCJj'
perc_NEUT_codes <- '^42b0|^XE2mP'
perc_LYMPH_codes <- '^42b1|^XE2mQ'
perc_MONO_codes <- '^42b2|^XE2mR'

ALP_codes <- '^44CU|^44F[.123Z]|^XaIRj|^XE2px'
ALT_codes <- '^44G[.3AB]|^X771f|^XaIRi|^XaLJx'
AST_codes <- '^44H[5BC]|^X771i|^XaES6|^XE25Q'
GGT_codes <- '^44G[479]|^XaES[34]'

sodium_codes <- '^44h[16]|^44I5|^4Q43|^X771T|^XaDva|^XaIRf|^XE2q0'
potassium_codes <- '^44h[08]|^44I4|^4Q42|^X771S|^XaDvZ|^XaIRl|^XE2pz'
calcium_codes <- '^44h[479]|^44I[8C]|^4Q721|^X771W|^Xabp[kr]|^XaDvd|^XaIR[kn]|^XaZyY|^XE2q3'
chloride_codes <- '^44h2|^44i1|^44I6|^XaDvb|^XaItN|^XE2q1'
phosphate_codes <- '^44h5|^44i2|^44I9|^XaDve|^XaItO|^XE2q4'
bicarbonate_codes <- '^44h3|^44i0|^44I7|^XaDvc|^XaItM|^XE2q2'
ferritin_codes <- '^42d4|^42R4|^XaItW|^XE24r'
iron_codes <- '^42d2|^42R[.1237Z]|^4Q74|^X76tH|^X7733|^XaIRe|^XE24q'
# Li_codes <- 'XE25g' # sample size is too small!

total_protein_codes <- '^44M[.1237AVZ]|^4Q3\\.|^XE2e[9C]' # omit XaJmG and XE25U
blood_globulin_codes <- '^44M5\\.|^44MP|^XaItX|^XE2eB'
blood_albumin_codes <- '^44M[4I]|^XaIRc|^XE2eA'

CRP_codes <- '^44C[CS]|^XaINL|^XE2dy'
PSA_codes <- '^43Z[2F]|^X80QD|^XabAM|^XaPqN|^XE25C'
CK_codes <- '^44H[4EG]|^XaES[5B]'
amylase_codes <- '^44C[4NT]|^XaIRh|^XM14L'
CA125_codes <- 'XaJNf'
# RF_codes <- '43F..' # hard to extract phenotype from this code.

Hct_codes <- '^425.|^X76t[bc]|^XE2Zq'
Hb_codes <- '^423[.123456789ABZ]|^Xa96v|^XE2m6|^XM1Vu'
MCH_codes <- '^428.|^XE2pb'
MCHC_codes <- '^429.'
MCV_codes <- '^42A[.12345Z]'
RDW_codes <- '^42Z7|^XE2mO'
ESR_codes <- '^42B6|^XE2m7'

MPV_codes <- '^42Z5'
PV_codes <- '^42B[.1Z]|^X76xg|^XE2pd'
INR_codes <- '^41C5|^42jr|^42QE[.01]|^9k21|^XaPJW'
# PT_codes <- '^42Q5' # usually use INR to replace PT: INR = PT/mean(normal PT), and hard to extract phenotype from this code.
APTT_codes <- '^42jG|^42Q6|^XS9U4'
# APTT_ratio_codes <- '^42Qu|^XaIU[24]' # sample size is too small!

urine_albumin_codes <- '^46N4|^46W[.1]|^XE2eI|^XE2bw'
urine_creatinine_codes <- '^46M7|^XE2qO'
UACR_codes <- '^44J7|^46T[CD]|^XE2n3'

blood_creatinine_codes <- '^44J3[.0123z]|^44J[CF]|^XE2q5|^XaETQ|^4Q40|^X771Q'
# eGFR_codes <- '^451[EFGKN]|^XacUK|^XaK8y|^XaMDA|^XaZpN|^XSFyN|^Y0a58'

bilirubin_codes <- '^44E[.1239CZ]|^XaERu|^XaETf|^XE2qu'
urea_codes <- '^44J[.1289AZ]|^X771P|^XaDvl|^XE25R|^XM0lt'
urate_codes <- '^44K|^4Q66|^XaDvn|^XE2e2|^XM0ls|^XM1Uq'
folate_codes <- '^42U[.1235EZ]|^X76tC'
B12_codes <- '^42T[.12Z]|^44L[ce]|^XaJ27|^XE2pf'
VD_codes <- '^44LA|^4QB4[67]|^XaY6m|^XE2e7'

TSH_codes <- '^442[AWX]|^X80Gc|^XaEL[VW]|^XE2wy'
FT3_codes <- '^442[45U]|^XaERq' # ignore 442f.
FT4_codes <- '^442[7cV]|^XaER[rs]'

# FSH_codes <- '^443[4hi]|^4Q23|^XaELZ|^XE25J|^XM0lx'
# LH_codes <- '^443[3ef]|^XaELa|^XE25I|^XM0lv'
# E2_codes <- '^446[57]|^XaEQo'
# T_codes <- '^447[3G]|^XaItQ|^XE2dr'
