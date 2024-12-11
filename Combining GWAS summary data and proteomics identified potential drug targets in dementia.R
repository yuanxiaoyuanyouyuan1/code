setwd("/data_analysis/zhaoyingjie/analysis")
library(data.table)
library(stringr)
library(TwoSampleMR)
library(dplyr)

# loading Exposure data ---------------------------------------------------

load('exposure.Rdata')
exp_dat<- exp_dat %>%
  dplyr::mutate(
    # Compute the F-statistic
    F_statistics = beta.exposure^2/se.exposure^2,
    # Compute the R^2
    pve = 2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure)
  )


# loading Outcome data ---------------------------------------------------

load('outcome.Rdata')



harmonised_plasma <- harmonise_data(plasma_exposure, Outcome)

harmonised_csf <- harmonise_data(CSF_exposure, Outcome)


# Main MR analysis ----------------------------------------------------------------


mr_csf <-TwoSampleMR::mr(harmonised_csf)
mr_plasma <- TwoSampleMR::mr(harmonised_plasma)

# Heterogeneity test --------------------------------------------

mr_heterogeneity_CSF <- mr_heterogeneity(harmonised_csf)
mr_heterogeneity_plasma <- mr_heterogeneity(harmonised_plasma)

mr_pleiotropy_test_CSF <- mr_pleiotropy_test(harmonised_csf)
mr_pleiotropy_test_plasma <- mr_pleiotropy_test(harmonised_plasma)


# MRPRESSO --------------------------------------------
library(MRPRESSO)

dat1 <- load('harmonized_data.Rdata')


mr_presso_result <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                                          SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                                          data = dat1, NbDistribution = 1000,  SignifThreshold = 0.05
)


# Bonferroni correction --------------------------------------------

load('Bonferroni_correction.Rdata')


disease <- unique(mr_csf$outcome)
i=1
c <- 0
table1_csf <- data.frame()
for(i in 1:length(disease)){
  a <- mr_csf[mr_csf$outcome%in%disease[i],]
  b <- length(  unique(a$exposure))
  h <- 0.05/b
  g <-  a%>%
    filter(pval < h,
           method %in% c("Wald ratio","Inverse variance weighted")) %>%
    left_join(pqtl_csf, by = "exposure")
  table1_csf <- rbind(table1_csf,g)
  c <- c+b

}


disease <- unique(mr_plasma$outcome)
i=1
c <- 0
table1_plasma <- data.frame()
for(i in 1:length(disease)){
  a <- mr_plasma[mr_plasma$outcome%in%disease[i],]
  b <- length(  unique(a$exposure))
  h <- 0.05/b
  g <-  a%>%
    filter(pval < h,
           method %in% c("Wald ratio","Inverse variance weighted")) %>%
    left_join(pqtl_csf, by = "exposure")
  table1_csf <- rbind(table1_csf,g)
  c <- c+b

}

# Summary of Sensitivity analysis --------------------------------

# Phenoscanner ---------------------------------------------------

mr_phenoscanner <- function(dat,
                            catalog = c("GWAS", "eQTL", "pQTL", "mQTL", "methQTL"),
                            pvalue = 5e-08,
                            proxies = "EUR",
                            r2 = 0.8,
                            build = 37)
{
  stopifnot("SNP" %in% names(dat))

  snpquery <- unique(dat$SNP)
  query_limit <- ceiling(0.01*length(snpquery))

  query_res <- list()
  for (i in 1:query_limit) {
    lo_i <- 100*i - 99
    up_i <- min(100*i,length(snpquery))
    snp_i <- snpquery[lo_i:up_i]
    print(paste0("Searching for ", snp_i))

    for (catalogue in catalog) {
      query_i <- try(
        phenoscanner::phenoscanner(
          snpquery = snp_i, catalogue = catalogue,
          proxies = proxies, pvalue = pvalue, r2 = r2, build = build)
      )
      if(!"try-error" %in% class(query_i))
      {
        if("results" %in% names(query_i))
        {
          query_i <- list(query_i$results)
          names(query_i) <- catalogue
          query_res <- append(query_res, query_i)
        }
      }
    }
  }

  return(query_res)
}



load('Phenoscanner.Rdata')

pleio_pheno <- mr_phenoscanner(Phenoscanner)

pleio_pheno_summary_raw <- pleio_pheno$GWAS %>%
  mutate(catalog = "GWAS") %>%
  bind_rows(pleio_pheno$pQTL%>%
              mutate(catalog = "pQTL")) %>%
  bind_rows(pleio_pheno$mQTL%>%
              mutate(catalog = "mQTL")) %>%
  filter(p != "NA",
         ancestry %in% c("Mixed", "European"))

g <- data.frame()
for (i in 1:nrow(pleio_pheno_summary_raw)){
  a <-pleio_pheno_summary_raw[i,]

  c <-b[b$snp %in% a$snp, ]
  f <- data.frame()
  if (nrow(c)==2){
    j=1
    for (j in 1:nrow(c)){
      e <- c[j,]
      d <- left_join(a,e,by="snp")
      f <- rbind(f,d)
    }

  }
  d <- left_join(a,c,by="snp")
  d <- rbind(d,f)
  g <- rbind(g,d)



}


h <- g %>%
  mutate(across(p, as.numeric)) %>%
  group_by(trait, pmid) %>%
  filter(p == min(p)) %>%
  ungroup() %>%
  dplyr::select(trait, Tissue,Protein = exposure, UniProt, snp, ref_a1, ref_a2, effect_allele.exposure:eaf.exposure, proxy, r2, rsid, a1, a2, pmid:p,catalog) %>%
  mutate(ref_a1f = ifelse(ref_a1 == effect_allele.exposure, eaf.exposure, 1- eaf.exposure),
         ref_a1f_reliability = ifelse((ref_a1 == effect_allele.exposure&ref_a2 == other_allele.exposure)|
                                        (ref_a2 == effect_allele.exposure&ref_a1 == other_allele.exposure),
                                      "high", "low")) %>%
  dplyr::select(trait, Tissue,Protein, UniProt, snp, effect_allele.exposure,ref_a1, ref_a2, ref_a1f,ref_a1f_reliability, proxy, r2, rsid, a1, a2, pmid:p,catalog) %>%
  rename(effector = effect_allele.exposure,
         SNP = rsid,
         effect_allele = a1,
         other_allele = a2) %>%
  unique() %>%
  rename(proxy_SNP = SNP,
         SNP = snp,
         pval = p) %>%
  mutate(effect_allele = effector,
         proxy_effect_allele = ifelse(effector == ref_a1, effect_allele, other_allele),
         beta = as.numeric(beta),
         se = as.numeric(se),
         beta = ifelse(effector == ref_a1, beta, -beta),
         proxy_SNP = ifelse(proxy==0,NA,proxy_SNP),
         r2 = ifelse(proxy==0,NA,r2),
         proxy_effect_allele = ifelse(proxy==0,NA,proxy_effect_allele)) %>%
  dplyr::select(Tissue:SNP,trait,catalog,effect_allele,proxy_SNP,r2,proxy_effect_allele,beta,se,pval,ancestry,pmid)

# Steiger filtering ----------------------------------------------

harmonised_csf %>%
  filter(SNP %in% table1_csf$SNP) %>%
  steiger_filtering()

harmonised_plasma %>%
  filter(SNP %in% table1_plasma$SNP) %>%
  steiger_filtering()


# Bidirectional MR ---------------------------------------------

load('bidirectional_exposure.Rdata')

load('bidirectional_outcome.Rdata')

harmonised_data <- harmonise_data(bidirectional_exposure, bidirectional_outcome)

mr_Bidirectional <- TwoSampleMR::mr(harmonised_data)


# Colocalization -----------------------------------------------
load('colocalization.Rdata')
library(coloc)

data <- list(dataset1 = pqtl, dataset2 = gwas)

p1 = 1e-04
p2 = 1e-04
p12 = 1e-05

res <- coloc::coloc.abf(dataset1 = data[[1]], dataset2 = data[[2]],
                        p1 = p1, p2 = p2, p12 = p12)



# MRlap for sample overlapping

library(MRlap)

overlapping = MRlap1(exposure = MRlap_exposure,
            exposure_name = "Protein",
            outcome = MRlap_Outcome,
            outcome_name = "Disease",
            ld = "~/eur_w_ld_chr",
            hm3 = "~/w_hm3.noMHC.snplist")



