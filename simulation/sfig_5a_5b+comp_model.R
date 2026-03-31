library(tidyr)
library(dplyr)
library(ggplot2)

# computational simulation of underlying biological process of resistance evolution
run_simulation = function(x, m, hA, hB, strategy){

  strategy = strategy
  
  # parameters
  x = x # overall compound exposure
  mA = m # compound efficacy on SS genotypes
  mB = m
  rA = 1 # fitness of RR genotypes for both alleles are 1 
  rB = 1
  hA = hA # dominance as determined by 2-day dose-response curves
  hB = hB
  cA=0 # both alleles are neutral, cost = 0
  cB=0
  hcA=0.5 # dominance of cost, set to 0.5, has no effect since cost is 0
  hcB=0.5
  
  #fitness
  wAABB = c(((1-x)+x*(1-mA*(1-rA)))*(1-cA)*(1-cB), 
            ((1-x)+x*(1-mB*(1-rB)))*(1-cA)*(1-cB),
            ((1-x)+x*(1-mA*(1-rA))*(1-mB*(1-rB)))*(1-cA)*(1-cB))
  wAaBB = c(((1-x)+x*(1-mA*(1-hA*rA)))*(1-hcA*cA)*(1-cB), 
            ((1-x)+x*(1-mB*(1-rB)))*(1-hcA*cA)*(1-cB), 
            ((1-x)+x*(1-mA*(1-hA*rA))*(1-mB*(1-rB)))*(1-hcA*cA)*(1-cB))
  waaBB = c(((1-x)+x*(1-mA))*(1-cB), 
            ((1-x)+x*(1-mB*(1-rB)))*(1-cB), 
            ((1-x)+x*(1-mA)*(1-mB*(1-rB)))*(1-cB))
  wAABb = c(((1-x)+x*(1-mA*(1-rA)))*(1-cA)*(1-hcB*cB), 
            ((1-x)+x*(1-mB*(1-hB*rB)))*(1-cA)*(1-hcB*cB), 
            ((1-x)+x*(1-mA*(1-rA))*(1-mB*(1-hB*rB)))*(1-cA)*(1-hcB*cB))
  wAaBb = c(((1-x)+x*(1-mA*(1-hA*rA)))*(1-hcA*cA)*(1-hcB*cB), 
            ((1-x)+x*(1-mB*(1-hB*rB)))*(1-hcA*cA)*(1-hcB*cB), 
            ((1-x)+x*(1-mA*(1-hA*rA))*(1-mB*(1-hB*rB)))*(1-hcA*cA)*(1-hcB*cB))
  waaBb = c(((1-x)+x*(1-mA))*(1-hcB*cB), 
            ((1-x)+x*(1-mB*(1-hB*rB)))*(1-hcB*cB), 
            ((1-x)+x*(1-mA)*(1-mB*(1-hB*rB)))*(1-hcB*cB))
  wAAbb = c(((1-x)+x*(1-mA*(1-rA)))*(1-cA), 
            ((1-x)+x*(1-mB))*(1-cA), 
            ((1-x)+x*(1-mA*(1-rA))*(1-mB))*(1-cA))
  wAabb = c((((1-x)+x*(1-mA*(1-hA*rA)))*(1-hcA*cA))*1, 
            ((1-x)+x*(1-mB))*(1-hcA*cA), 
            ((1-x)+x*(1-mA*(1-hA*rA))*(1-mB))*(1-hcA*cA))
  waabb = c(((1-x)+x*(1-mA)), 
            ((1-x)+x*(1-mB)), 
            ((1-x)+x*(1-mA)*(1-mB)))
  fitness = list(wAABB, wAaBB, waaBB, wAABb, wAaBb, waaBb, wAAbb, wAabb, waabb)

  # simulation setup
  t = 15
  reps = 500
  k = 3500 # population size 
  fA0 = 0.02 # initial allele frequency of resistance to insecticide A (hets)
  fB0 = 0.02 # initial allele frequency of resistance to insecticide B (hets)
  
  # setup of data to save 
  fA_reps = matrix(fA0/2,t+1,reps) # store for resistant allele frequency across reps
  fB_reps = matrix(fB0/2,t+1,reps) # store for resistant allele frequency across reps
  
  # simulation 
  for(i in 1:reps){
    for(j in 1:t){
      # starting values 
      if(j==1){ 
        fA = fA0 # resistant allele frequency 
        fB = fB0 # resistant allele frequency 
        n = k # population size 
        # linkage disequilibrium
        genotypes = matrix(0,t+1,9)
        colnames(genotypes) = c("fAABB","fAaBB","faaBB",
                                "fAABb","fAaBb","faaBb",
                                "fAAbb","fAabb","faabb")
        genotypes[1,] = c(0,0,0,0,0,fB0,0,fA0,1-fA0-fB0) # starting frequency is half the fA0 value due to startin in hets
      }
      
      # genotype frequencies
      fAABB = genotypes[j,1]
      fAaBB = genotypes[j,2]
      faaBB = genotypes[j,3]
      fAABb = genotypes[j,4]
      fAaBb = genotypes[j,5]
      faaBb = genotypes[j,6]
      fAAbb = genotypes[j,7]
      fAabb = genotypes[j,8]
      faabb = genotypes[j,9]
      
      # choose strategy
      if(strategy==1){ # A only
        a = 1
      } else if(strategy == 2){ # B only
        a = 2
      } else if(strategy == 3){ # Mixture
        a = 3
      } else if(strategy == 4){ # Rotation
        if(j%%2 == 1){ # odd generations
          a = 1 # use A
        } else if(j%%2 == 0){ # even generations
          a = 2 # use B
        }
      }
      # gamete frequencies 
      fAB = wAABB[a]*fAABB + 0.5*wAaBB[a]*fAaBB + 0.5*wAABb[a]*fAABb + 0.25*wAaBb[a]*fAaBb
      faB = 0.5*wAaBB[a]*fAaBB + waaBB[a]*faaBB + 0.25*wAaBb[a]*fAaBb + 0.5*waaBb[a]*faaBb
      fAb = 0.5*wAABb[a]*fAABb + 0.25*wAaBb[a]*fAaBb + wAAbb[a]*fAAbb + 0.5*wAabb[a]*fAabb
      fab = 0.25*wAaBb[a]*fAaBb + 0.5*waaBb[a]*faaBb + 0.5*wAabb[a]*fAabb + waabb[a]*faabb
      # proportional gamete frequencies, sums to one 
      wm = fAB+faB+fAb+fab
      pAB = fAB/wm
      paB = faB/wm
      pAb = fAb/wm
      pab = fab/wm
      # check: pAB + paB + pAb + pab
      # zygote frequencies 
      pAABB = pAB*pAB
      pAaBB = 2*pAB*paB
      paaBB = paB*paB
      pAABb = 2*pAB*pAb 
      pAaBb = 2*pAB*pab + 2*pAb*paB
      paaBb = 2*paB*pab 
      pAAbb = pAb*pAb
      pAabb = 2*pAb*pab
      paabb = pab*pab
      
      # stochastic sampling
      pset = c(pAABB,pAaBB,paaBB,pAABb,pAaBb,paaBb,pAAbb,pAabb,paabb)
      qset = rmultinom(1,k,prob=pset) # sample from multinomial distribution 
      pAABB = qset[1]/k
      pAaBB = qset[2]/k
      paaBB = qset[3]/k
      pAABb = qset[4]/k
      pAaBb = qset[5]/k
      paaBb = qset[6]/k
      pAAbb = qset[7]/k
      pAabb = qset[8]/k
      paabb = qset[9]/k
      
      # proportional allele frequencies 
      pA = pAABB + 0.5*pAaBB + pAABb + 0.5*pAaBb + pAAbb + 0.5*pAabb
      pa = 0.5*pAaBB + paaBB + 0.5*pAaBb + paaBb + 0.5*pAabb + paabb
      # check: pA + pa
      pB = pAABB + pAaBB + paaBB + 0.5*pAABb + 0.5*pAaBb + 0.5*paaBb
      pb = 0.5*pAABb + 0.5*pAaBb + 0.5*paaBb + pAAbb + pAabb + paabb
      # check: pB + pb
      # data to save 
      fA_reps[1+j,i] = pA # resistance allele frequency gets saved to the main matrix
      fB_reps[1+j,i] = pB
      
      genotypes[1+j,1] = pAABB
      genotypes[1+j,2] = pAaBB
      genotypes[1+j,3] = paaBB
      genotypes[1+j,4] = pAABb
      genotypes[1+j,5] = pAaBb
      genotypes[1+j,6] = paaBb
      genotypes[1+j,7] = pAAbb
      genotypes[1+j,8] = pAabb
      genotypes[1+j,9] = paabb
      
      # update values
      fA = pA
      fB = pB
    }
  }
  
  #data to save
  generations = seq(0, t)
  
  df_A = as.data.frame(fA_reps)
  df_A$Generation = generations
  df_A_long = df_A %>%
    pivot_longer(-Generation, names_to = "Replicate", values_to = "Frequency")
  
  df_B = as.data.frame(fB_reps)
  df_B$Generation = generations
  df_B_long = df_B %>%
    pivot_longer(-Generation, names_to = "Replicate", values_to = "Frequency")
  
  out = list(A_freq = df_A_long, B_freq = df_B_long, fitness = fitness)
  return(out)
}

# dominance calculation from empirical data

# carbendazim
# at 7.7 ug/ml, SS fitness = 0.101, RS fitness = 0.635, RR fitness = 1
# at 6.1 ug/ml, SS fitness = 0.308, RS fitness = 0.791, RR fitness = 1
# spirotetramat
# at 46 ug/ml, SS fitness = 0.0993, RS fitness = 0.620, RR fitness = 1
# at 33.5 ug/ml, SS fitness = 0.302, RS fitness = 0.701, RR fitness = 1

carb_h_high = (0.635-0.101)/(1-0.101)
carb_h_low = (0.791-0.308)/(1-0.308)
spir_h_high = (0.62-0.0993)/(1-0.0993)
spir_h_low = (0.701-0.302)/(1-0.302)

# run simulations at experimentally tested conditions
data1_high = run_simulation(x=0.9,m=0.9,hA=carb_h_high,hB=spir_h_high,strategy=3)
data1_low = run_simulation(x=0.7,m=0.9,hA=carb_h_high,hB=spir_h_high,strategy=3)

data2_high = run_simulation(x=0.9,m=0.7,hA=carb_h_low,hB=spir_h_low,strategy=3)
data2_low = run_simulation(x=0.7,m=0.7,hA=carb_h_low,hB=spir_h_low,strategy=3)

data3_high = run_simulation(x=0.9,m=0.9,hA=carb_h_high,hB=spir_h_high,strategy=4)
data3_low = run_simulation(x=0.7,m=0.9,hA=carb_h_high,hB=spir_h_high,strategy=4)

data4_high = run_simulation(x=0.9,m=0.7,hA=carb_h_low,hB=spir_h_low,strategy=4)
data4_low = run_simulation(x=0.7,m=0.7,hA=carb_h_low,hB=spir_h_low,strategy=4)

# two versions of function to summarise allele A and B data
A_stats = function(data) {
  data$A_freq %>%
    group_by(Generation) %>%
    summarise(
      mean = mean(Frequency, na.rm = TRUE),
      p2.5 = quantile(Frequency, 0.025, na.rm = TRUE),
      p97.5 = quantile(Frequency, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}
A_stats = function(data) {
  data$B_freq %>%
    group_by(Generation) %>%
    summarise(
      mean = mean(Frequency, na.rm = TRUE),
      p2.5 = quantile(Frequency, 0.025, na.rm = TRUE),
      p97.5 = quantile(Frequency, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

data1_high_A = A_stats(data1_high)
data1_low_A = A_stats(data1_low)
data2_high_A = A_stats(data2_high)
data2_low_A = A_stats(data2_low)
data3_high_A = A_stats(data3_high)
data3_low_A = A_stats(data3_low)
data4_high_A = A_stats(data4_high)
data4_low_A = A_stats(data4_low)

data1_high_A$condition = "A"
data1_low_A$condition = "F"
data2_high_A$condition = "D"
data2_low_A$condition = "H"
data3_high_A$condition = "B"
data3_low_A$condition = "G"
data4_high_A$condition = "E"
data4_low_A$condition = "I"

data = rbind(data1_high_A, data1_low_A,
             data2_high_A, data2_low_A,
             data3_high_A, data3_low_A,
             data4_high_A, data4_low_A)

data = data %>%
  mutate(
    regimen = if_else(condition %in% c("A", "D", "F", "H"),
                      "Mixture",
                      "Rotation")
  ) 
data = data %>%
  mutate(
    efficacy = if_else(condition %in% c("A", "B", "F", "G"),
                       "90% efficacy",
                       "70% efficacy")
  )
data = data %>%
  mutate(
    exposure = if_else(condition %in% c("A", "B", "D", "E"),
                       "90% exposure",
                       "70% exposure")
  )

data$regimen = factor(data$regimen, levels=c("Mixture", "Rotation"))
data$efficacy = factor(data$efficacy, levels=c("90% efficacy", "70% efficacy"))
data$exposure = factor(data$exposure, levels=c("90% exposure", "70% exposure"))

ggplot(data, aes(Generation, mean)) +
  geom_line(linewidth = 0.6, aes(color = exposure)) +
  geom_ribbon(aes(x = Generation, ymin = p2.5, ymax = p97.5, group = exposure),fill = "grey50",alpha = 0.1)+
  facet_grid(regimen ~ efficacy) +
  theme_bw(base_size=12) +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  labs(x = "Generation", y = "Pred. spirotetramat resistance allele frequency", color = "Exposure") +
  scale_color_manual(values = c("90% exposure" = "#6F6F6F", "70% exposure" = "#C65858")) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    legend.position = "bottom"
  )

ggsave(
  filename = "sfig_5a.pdf",
  width = 170, height = 135,
  units = "mm",
  dpi = 1200,
  device = cairo_pdf
)






