######################################################
# Estimate SNP heritability of breast cancer and SCZ #
######################################################

cd /users/zhe/Adv-genetic-stats-assignment-2024/assignment_1

module load anaconda
module load ldsc

wget --recursive --no-parent https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/

# Munge data
munge_sumstats.py \
    --sumstats INPUT/data_BC17.txt \
    --out INPUT/BC17

munge_sumstats.py \
    --sumstats INPUT/data_SCZ.txt \
    --ignore beta \
    --out INPUT/SCZ

# LDSC
ldsc.py \
    --h2 INPUT/BC17.sumstats.gz \
    --ref-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --w-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --out OUTPUT/BC17_bip
cat OUTPUT/BC17_bip.log

ldsc.py \
    --h2 INPUT/SCZ.sumstats.gz \
    --ref-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --w-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --out OUTPUT/SCZ_bip
cat OUTPUT/SCZ_bip.log

# LDSC (liability adjustment)
ldsc.py \
    --h2 INPUT/BC17.sumstats.gz \
    --ref-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --w-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --out OUTPUT/BC17_bip_lia \
    --samp-prev 0.5267109 \
    --pop-prev 0.00059
cat OUTPUT/BC17_bip_lia.log

ldsc.py \
    --h2 INPUT/SCZ.sumstats.gz \
    --ref-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --w-ld-chr ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ \
    --out OUTPUT/SCZ_bip_lia \
    --samp-prev 0.4309786 \
    --pop-prev 0.0032
cat OUTPUT/SCZ_bip_lia.log

# Also can transform in R
module load conda_R
R

K = 0.00059
P = 0.5267109
h02 = 0.5541
c = qnorm(1-K)
z = dnorm(c)
h02*(K*(1-K)/z^2)*(K*(1-K)/P/(1-P))

K = 0.0032
P = 0.4309876
h02 = 1.426
c = qnorm(1-K)
z = dnorm(c)
h02*(K*(1-K)/z^2)*(K*(1-K)/P/(1-P))

quit(save = "no")

# AUC
R

pnorm(sqrt(0.5541)/sqrt(2))
pnorm(sqrt(1.426)/sqrt(2))

pnorm(sqrt(0.5541*0.5)/sqrt(2))
pnorm(sqrt(1.426*0.5)/sqrt(2))

quit(save = "no")
