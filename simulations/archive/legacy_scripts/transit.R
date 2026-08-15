a <- function(rho, p, q) {
  rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
}
n <- 1000
p <- 0.2 
q <- 0.6
rho <- 0.25
p_cases_exp <-  rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
prob <- c(`p1`=cont,
          `p2`=1-q-cont,
          `p3`=1-p-cont,
          `p4`=cont+p+q-1)
n.sim <- 1
u <- sample.int(4, n, replace=TRUE, prob=prob)
y <- floor((u-1)/2)
x <- 1 - u %% 2 # 1 - u %% 2

tab = table(x,y)
n_cases_exp = tab[2,1]
n_controls_exp = tab[2,2]
n_cases_nexp = tab[1,1]
n_controls_nexp = tab[1,2]

n_cases = n_cases_exp+n_cases_nexp
n_controls = n_controls_exp+n_controls_nexp
n_exp = n_cases_exp + n_controls_exp
n_nexp = n_cases_nexp + n_controls_nexp

r_raw_data = cor.test(~x+y)$estimate
n_exp;n_nexp
n_cases;n_controls
tab
mean(x); mean(y)
n_exp/(n_nexp+n_exp)
n_cases/(n_nexp+n_exp)
cor.test(~x+y)$estimate
