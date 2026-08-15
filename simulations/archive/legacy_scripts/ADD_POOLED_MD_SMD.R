library(MonteCarlo); library(tidyverse); library(ggplot2)
d <- c(70, 85, 100, 115, 130)[1]
het <- c(0, 10, 20)[1]
n_study <- c(5, 10, 20)[1]
n_sample <- c(20, 50, 100)[1]
d1 <- c(70, 85, 100, 115, 130)
het1 <- c(0, 10, 20, 30)
n_study1 <- c(5, 10, 20)
n_sample1 <- c(20, 50, 100)

dat = data.frame(d = NA, het = NA,
                 n_study = NA, n_sample = NA,
                 I2 = NA,
                 smd_true = NA,
                 smd_est = NA,
                 pval_true = NA,
                 pval_est = NA)

for (d in d1) {
  for (het in het1) {
    for (n_study in n_study1) {
      for (n_sample in n_sample1) {
        for (i in 1:10) {
          dat2 = smd_true = smd_est = pval_true = pval_est = NULL

          mean_nexp = runif(n_study, 100 - het, 100 + het)
          sd_exp = sd_nexp = rep(15, n_study)
          mean_exp = rnorm(n_study, d, 1)
          n_exp = n_nexp = rep(n_sample/2, n_study)

          res_md = metafor::rma.uni(
            m1i = mean_exp, m2i = mean_nexp,
            sd1i = sd_exp, sd2i = sd_nexp,
            n1i = n_exp, n2i = n_nexp,
            measure = "MD")

          res_smd = metafor::rma.uni(
            m1i = mean_exp, m2i = mean_nexp,
            sd1i = sd_exp, sd2i = sd_nexp,
            n1i = n_exp, n2i = n_nexp,
            measure = "SMD")

          res_md_sd = res_md$se / sqrt(1/sum(n_exp) + 1/sum(n_nexp))
          smd_est = as.numeric(res_md$beta) / res_md_sd
          smd_true = as.numeric(res_smd$beta)

          pval_est = res_md$pval
          pval_true = res_smd$pval
          dat2 = data.frame(d = d, het = het,
                            n_study = n_study,
                            n_sample = n_sample,
                            I2 = res_md$I2,
                            smd_true,
                            smd_est,
                            pval_true,
                            pval_est)
          dat = rbind(dat, dat2)
        }
      }
    }
  }
}

dat$I2_cat <- cut(dat$I2,
     breaks = c(-Inf,25,50, 75,Inf),
     labels = c("I² = 0-25%", "I² = 26-50%",
                "I² = 51-75%", "I² = 76-100%"))
dat = dat[2:nrow(dat), ]

dat$smd_bias = dat$smd_true - dat$smd_est
dat$pval_bias = dat$pval_true - dat$pval_est

res = dat %>%
  group_by(d, I2_cat) %>%
  summarise(smd_bias = mean(smd_bias),
            pval_bias = mean(pval_bias))

res$d[res$d == "70"] <-  "SMD = -2"
res$d[res$d == "85"] <-  "SMD = -1"
res$d[res$d == "100"] <- "SMD = 0"
res$d[res$d == "115"] <- "SMD = 1"
res$d[res$d == "130"] <- "SMD = 2"

res$d <- factor(res$d,
                levels = c("SMD = -2",
                           "SMD = -1",
                           "SMD = 0",
                           "SMD = 1",
                           "SMD = 2"))

ggplot(res, aes(x = d, y = smd_bias)) +
  geom_bar(stat = "identity",
           aes(fill = d), alpha = 0.5) +
  facet_wrap(~ I2_cat) +
  theme_bw()


