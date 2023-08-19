VCrelchange_entomology_GSA <- function(HBI_mean, HBI_samples, parous_mean, parous_samples,
                                   sac_mean, sac_samples, tau_mean, tau_samples,
                                   N_b){
  theta_d = 8/24
  P_B_b = 0.95
  P_C_b = 0.95
  P_D_b = 0.99
  P_E_b = 0.88
  theta_s = 10
  eta = 1
  #baseline
  chi = HBI_mean
  Pf = parous_mean
  A0 = sac_mean
  tau = tau_mean
  PA = 1 - A0
  PA1 = A0*Pf*chi*eta*P_B_b*P_C_b/(P_B_b*P_C_b*(chi*P_D_b*P_E_b + (1-chi)*P_D_b*P_E_b)*eta*P_B_b*P_C_b)
  PAh = A0*Pf*(1-chi)/(eta*P_B_b*P_C_b*(chi*P_D_b*P_E_b +(1-chi)*P_D_b*P_E_b))
  PA2 = eta*PAh
  alpha_H_b = (1/N_b)*(PA1/(1-PA))*(-log(PA)/theta_d)
  alpha_A_b = (1/N_b)*(PA2/(1-PA))*(-log(PA)/theta_d)
  mu_A_b = ((1-(PA + PA1 +PA2))/(1-PA))*(-log(PA)/theta_d)
  K_plus = floor(theta_s/tau) - 1
  #PA = exp(-(alpha_H_b*N_b + mu_A_b)*theta_d)
  #PAi = (1 - PA)*((alpha_H_b*N_b)/(alpha_H_b*N_b + mu_A_b))
  #Pdf = PAi*P_B_b*P_C_b*P_D_b*P_E_b
  MalariaHost = c(1,0)
  Hosts = which(MalariaHost == 1)
  N = c(N_b, N_b)
  alpha_H = c(alpha_H_b, alpha_A_b)
  P_C = c(P_C_b, P_C_b)
  PA = exp(-(sum(alpha_H*N) + mu_A_b)*theta_d)
  PAi = (1 - PA)*((alpha_H*N)/(sum(alpha_H*N) + mu_A_b))
  #PA + sum(PAi) + (1 - PA)*(mu_A_b/(sum(alpha_H*N) + mu_A_b))
  Pdf = sum(PAi*P_B_b*P_C*P_D_b*P_E_b)
  #####
  frac1 = sum(PAi[Hosts]*P_B_b)*sum(PAi[Hosts]*P_B_b*P_C[Hosts]*P_D_b*P_E_b)/((1-PA-Pdf)^2)
  frac2 = 0
  for(h in 0:K_plus){
    frac2 = frac2 + choose(theta_s - (h+1)*tau + h, h)*(PA^(theta_s - (h+1)*tau))*Pdf^h
  }
  frac3 = 0
  for(l in 1:(tau-1)){
    for(h in 0:(floor((theta_s + l)/tau) - 2)){
      frac3 = frac3 + choose(theta_s + l - (h+2)*tau + h, h)*(PA^(theta_s + l -(h+2)*tau))*Pdf^(h+1)
    }
  }
  baseline = frac1*(frac2 + frac3)/N_b
  ## relative change
    chi = HBI_samples
    Pf = parous_samples
    A0 = sac_samples
    tau = tau_samples
    PA = 1 - A0
    PA1 = A0*Pf*chi*eta*P_B_b*P_C_b/(P_B_b*P_C_b*(chi*P_D_b*P_E_b + (1-chi)*P_D_b*P_E_b)*eta*P_B_b*P_C_b)
    PAh = A0*Pf*(1-chi)/(eta*P_B_b*P_C_b*(chi*P_D_b*P_E_b +(1-chi)*P_D_b*P_E_b))
    PA2 = eta*PAh
    alpha_H_b = (1/N_b)*(PA1/(1-PA))*(-log(PA)/theta_d)
    alpha_A_b = (1/N_b)*(PA2/(1-PA))*(-log(PA)/theta_d)
    mu_A_b = ((1-(PA + PA1 +PA2))/(1-PA))*(-log(PA)/theta_d)
    K_plus = floor(theta_s/tau) - 1
    #PA = exp(-(alpha_H_b*N_b + mu_A_b)*theta_d)
    #PAi = (1 - PA)*((alpha_H_b*N_b)/(alpha_H_b*N_b + mu_A_b))
    #Pdf = PAi*P_B_b*P_C_b*P_D_b*P_E_b
    MalariaHost = c(1,0)
    Hosts = which(MalariaHost == 1)
    N = c(N_b, N_b)
    alpha_H = c(alpha_H_b, alpha_A_b)
    P_C = c(P_C_b, P_C_b)
    PA = exp(-(sum(alpha_H*N) + mu_A_b)*theta_d)
    PAi = (1 - PA)*((alpha_H*N)/(sum(alpha_H*N) + mu_A_b))
    #PA + sum(PAi) + (1 - PA)*(mu_A_b/(sum(alpha_H*N) + mu_A_b))
    Pdf = sum(PAi*P_B_b*P_C*P_D_b*P_E_b)
    #####
    # sherry check h in frac 3, should it be tau?
    frac1 = sum(PAi[Hosts]*P_B_b)*sum(PAi[Hosts]*P_B_b*P_C[Hosts]*P_D_b*P_E_b)/((1-PA-Pdf)^2)
    frac2 = 0
    for(h in 0:K_plus){
      frac2 = frac2 + choose(theta_s - (h+1)*tau + h, h)*(PA^(theta_s - (h+1)*tau))*Pdf^h
    }
    frac3 = 0
    for(l in 1:(tau-1)){
      for(h in 0:(floor((theta_s + l)/tau) - 2)){
        frac3 = frac3 + choose(theta_s + l - (h+2)*tau + h, h)*(PA^(theta_s + l -(h+2)*tau))*Pdf^(h+1)
      }
    }
    other = frac1*(frac2 + frac3)/N_b
    RelativeChange = (other-baseline)/baseline
  return(RelativeChange)
}
