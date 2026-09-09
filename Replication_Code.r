################################################################################
# REPLICATION CODE FOR: "Optimal longevity of a dynasty"
# AUTHORS: Satoshi Nakano and Kazuhiko Nishimura
# DESCRIPTION: This script generates the theoretical figures (Fig 2 - Fig 6).
################################################################################

# --- Initial Setup & Package Loading ---
# Uncomment the following lines to install required packages if not already installed
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("DescTools")
# install.packages("gglorenz")

library(ggplot2)
library(dplyr)
library(DescTools)
library(gglorenz)

################################################################################
# FIGURE 2: AK Setting (First Order Condition & Population Value)
################################################################################
cat("Generating Figure 2...\n")
beta = 0.992
alpha1 = 1.012
alpha2 = 1.010
alpha3 = 1/beta
alpha4 = 1.005
k0 = 150
H = 600

# Population Value Data
V0_Ak1 <- numeric(); V0_Ak2 <- numeric(); V0_Ak3 <- numeric(); V0_Ak4 <- numeric()
for (h in 1:(H+1)) {
  V0_Ak1[h] = log(alpha1*beta)*(beta-((1-beta)*(h-1)+1)*beta^(h))/((1-beta)^2)+((1-beta^(h))/(1-beta))*log(alpha1*k0*(1-beta)/(1-beta^(h)))
  V0_Ak2[h] = log(alpha2*beta)*(beta-((1-beta)*(h-1)+1)*beta^(h))/((1-beta)^2)+((1-beta^(h))/(1-beta))*log(alpha2*k0*(1-beta)/(1-beta^(h)))
  V0_Ak3[h] = log(alpha3*beta)*(beta-((1-beta)*(h-1)+1)*beta^(h))/((1-beta)^2)+((1-beta^(h))/(1-beta))*log(alpha3*k0*(1-beta)/(1-beta^(h)))
  V0_Ak4[h] = log(alpha4*beta)*(beta-((1-beta)*(h-1)+1)*beta^(h))/((1-beta)^2)+((1-beta^(h))/(1-beta))*log(alpha4*k0*(1-beta)/(1-beta^(h)))
}
TH <- seq(0,H,length=H+1)
df_Ak <- data.frame(V0_Ak1, V0_Ak2, V0_Ak3, V0_Ak4, TH)

ggplot() + 
  geom_line(data=df_Ak, aes(x=TH, y=V0_Ak1), linetype="solid") +
  geom_line(data=df_Ak, aes(x=TH, y=V0_Ak2), linetype="solid") +
  geom_line(data=df_Ak, aes(x=TH, y=V0_Ak3), linetype="solid") +
  geom_line(data=df_Ak, aes(x=TH, y=V0_Ak4), linetype="solid") +
  theme(aspect.ratio = 1) +  
  xlab("Planning Horizon") + ylab("Population Value")
ggsave(file = "Figure_2_PV.eps", width = 6.56, height = 4.68)

# FOC Data
FOC1 <- numeric(); FOC2 <- numeric(); FOC3 <- numeric(); FOC4 <- numeric(); FOCC <- numeric()
PHI1 = -1 + (1/(1-beta) + (1/log(beta)))*log(alpha1*beta) +log(alpha1*k0*(1-beta))
PHI2 = -1 + (1/(1-beta) + (1/log(beta)))*log(alpha2*beta) +log(alpha2*k0*(1-beta))
PHI3 = -1 + (1/(1-beta) + (1/log(beta)))*log(alpha3*beta) +log(alpha3*k0*(1-beta))
PHI4 = -1 + (1/(1-beta) + (1/log(beta)))*log(alpha4*beta) +log(alpha4*k0*(1-beta))

for (n in 1:(H+1)) {
  FOC1[n] = ((n-1)*log(alpha1*beta) + PHI1)
  FOC2[n] = ((n-1)*log(alpha2*beta) + PHI2)
  FOC3[n] = ((n-1)*log(alpha3*beta) + PHI3)
  FOC4[n] = ((n-1)*log(alpha4*beta) + PHI4)
  FOCC[n] = (log(1-beta^n))
}
df_FOC <- data.frame(FOC1, FOC2, FOC3, FOC4, FOCC, TH)

ggplot() + 
  geom_line(data=df_FOC, aes(x=TH, y=FOC1), linetype="solid") +
  geom_line(data=df_FOC, aes(x=TH, y=FOC2), linetype="solid") +
  geom_line(data=df_FOC, aes(x=TH, y=FOC3), linetype="solid") +
  geom_line(data=df_FOC, aes(x=TH, y=FOC4), linetype="solid") +
  geom_line(data=df_FOC, aes(x=TH, y=FOCC), linetype="dashed") +
  theme(aspect.ratio = 1) +
  xlab("Planning Horizon") + ylab("First Order Condition")
ggsave(file = "Figure_2_FOC.eps", width = 6.56, height = 4.68)


################################################################################
# FIGURE 3: AK Setting (Trajectories of Undiscounted Contribution & Capital)
################################################################################
cat("Generating Figure 3...\n")
alpha = 1.012 
lnC_200 <- numeric(); b_lnC_200 <- numeric(); k_200 <- numeric()
lnC_400 <- numeric(); b_lnC_400 <- numeric(); k_400 <- numeric()
lnC_600 <- numeric(); b_lnC_600 <- numeric(); k_600 <- numeric()

n_200 = 200
for (t in 1:(n_200+1)) {
  lnC_200[t] = log(alpha*k0*(alpha*beta)^(t-1)) -log((1-beta^(n_200+1))/(1-beta))
  b_lnC_200[t] = (beta^t)*lnC_200[t] 
  k_200[t] = (1/alpha)*exp(lnC_200[t])*(1-beta^(n_200+1-(t-1)))/(1-beta)
}
n_400 = 400
for (t in 1:(n_400+1)) {
  lnC_400[t] = log(alpha*k0*(alpha*beta)^(t-1)) -log((1-beta^(n_400+1))/(1-beta))
  b_lnC_400[t] = (beta^t)*lnC_400[t] 
  k_400[t] = (1/alpha)*exp(lnC_400[t])*(1-beta^(n_400+1-(t-1)))/(1-beta)
}
n_600 = 600
for (t in 1:(n_600+1)) {
  lnC_600[t] = log(alpha*k0*(alpha*beta)^(t-1)) -log((1-beta^(n_600+1))/(1-beta))
  b_lnC_600[t] = (beta^t)*lnC_600[t] 
  k_600[t] = (1/alpha)*exp(lnC_600[t])*(1-beta^(n_600+1-(t-1)))/(1-beta)
}

th_200 <- seq(0, n_200, length=n_200+1)
th_400 <- seq(0, n_400, length=n_400+1)
th_600 <- seq(0, n_600, length=n_600+1)
df_lnC_200 <- data.frame(lnC_200, k_200, th_200)
df_lnC_400 <- data.frame(lnC_400, k_400, th_400)
df_lnC_600 <- data.frame(lnC_600, k_600, th_600)

ggplot() + 
  geom_line(data=df_lnC_200, aes(x=th_200, y=lnC_200), linetype="solid") +
  geom_line(data=df_lnC_400, aes(x=th_400, y=lnC_400), linetype="solid") +
  geom_line(data=df_lnC_600, aes(x=th_600, y=lnC_600), linetype="solid") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme(aspect.ratio = 1) +  
  xlab("Generation") + ylab("Contribution (Undiscounted)")
ggsave(file = "Figure_3_Contribution.eps", width = 6.56, height = 4.68)

ggplot() + 
  geom_line(data=df_lnC_200, aes(x=th_200, y=k_200), linetype="solid") +
  geom_line(data=df_lnC_400, aes(x=th_400, y=k_400), linetype="solid") +
  geom_line(data=df_lnC_600, aes(x=th_600, y=k_600), linetype="solid") +
  theme(aspect.ratio = 1) +  
  xlab("Generation") + ylab("Capital Intensity")
ggsave(file = "Figure_3_Capital.eps", width = 6.56, height = 4.68)


################################################################################
# FIGURE 4: ZD Setting (Population Value & Trajectories)
################################################################################
cat("Generating Figure 4...\n")
beta = 0.99999999999
theta = 1 
alpha = 1.0001
H = 600

# Population Value
V <- numeric(); pvlnC2 <- numeric()
for (h in 1:H) {
  for (t in 1:h) {
    lnB2 <- numeric()
    for (i in 1:t) {
      lnB2[i] = (theta^(i))*log(theta*alpha*beta*(1-(beta*theta)^(h-t+i))/(1-(beta*theta)^(h-t+i+1)))
    }
    lnC2_val = log(alpha) + (theta^(t+1))*log(k0) - log((1-(beta*theta)^(h-t+1))/(1-beta*theta)) + sum(lnB2[1:t])
    pvlnC2[t] = (beta^(t))*lnC2_val
  }
  V[h] = log((alpha*(k0)^theta)/((1-(beta*theta)^(h+1))/(1-beta*theta))) + sum(pvlnC2[1:h])
}
V0 <- c(log((alpha*(k0)^theta)), V)
TH <- seq(0, H, length=H+1)
df_ZD <- data.frame(V0, TH)

ggplot(df_ZD, aes(x = TH, y = V0)) + 
  geom_line(linetype="solid") +
  theme(aspect.ratio = 1) +  
  xlab("Planning Horizon") + ylab("Population Value")
ggsave(file = "Figure_4_PV.eps", width = 6.56, height = 4.68)


################################################################################
# FIGURE 5: Lorenz Curves
################################################################################
cat("Generating Figure 5...\n")
theta = 0.955392
beta = 1
alpha = 1.2
k0 = 150

generate_lorenz_data <- function(H_val, id_label) {
  Co <- numeric()
  for (t in 1:H_val) {
    lnB2 <- numeric()
    for (i in 1:t) {
      lnB2[i] = (theta^(i))*log(theta*alpha*beta*(1-(beta*theta)^(H_val-t+i))/(1-(beta*theta)^(H_val-t+i+1)))
    }
    lnC2_val = log(alpha) + (theta^(t+1))*log(k0) - log((1-(beta*theta)^(H_val-t+1))/(1-beta*theta)) + sum(lnB2[1:t])
    Co[t] = exp(lnC2_val)
  }
  data.frame(id = rep(id_label, H_val), Co = Co)
}

df_a <- generate_lorenz_data(200, 'a')
df_b <- generate_lorenz_data(400, 'b')
df_c <- generate_lorenz_data(600, 'c')
df_abc <- dplyr::bind_rows(df_a, df_b, df_c)

ggplot(df_abc, aes(x = Co, linetype = id)) +
  stat_lorenz() +
  theme(legend.position = 'none') +
  theme(aspect.ratio = 1) +
  geom_abline(color="grey") +
  xlab("Cumulative Population") + ylab("Cumulative Share")
ggsave(file = "Figure_5_Lorenz.eps", width = 6.56, height = 4.68)


################################################################################
# FIGURE 6: Gini Index
################################################################################
cat("Generating Figure 6...\n")
theta = 0.99
beta = 0.9
alpha = 1  
k0 = 150
H = 600

G <- numeric()
for (h in 1:H) {
  Co <- numeric()
  for (t in 1:h) {
    lnB2 <- numeric()
    for (i in 1:t) {
      lnB2[i] = (theta^(i))*log(theta*alpha*beta*(1-(beta*theta)^(h-t+i))/(1-(beta*theta)^(h-t+i+1)))
    }
    lnC2_val = log(alpha) + (theta^(t+1))*log(k0) - log((1-(beta*theta)^(h-t+1))/(1-beta*theta)) + sum(lnB2[1:t])
    Co[t] = exp(lnC2_val)
  }
  G[h] = Gini(Co)
}
G_a <- G  # Solid line data
G_b <- G  # Placeholder for dashed line data (run with different parameters if desired)

TH_Gini <- seq(1, H, length=H)
df_Gini <- data.frame(G_a, G_b, TH_Gini)

ggplot(df_Gini) + 
  geom_line(aes(x=TH_Gini, y=G_a), linetype="solid") +
  geom_line(aes(x=TH_Gini, y=G_b), linetype="dashed") +
  theme(aspect.ratio = 1) +  
  xlab("Planning Horizon") + ylab("Gini Index")
ggsave(file = "Figure_6_Gini.eps", width = 6.56, height = 4.68)

cat("All figures generated successfully!\n")