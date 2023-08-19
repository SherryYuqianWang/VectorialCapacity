VCrelchange_intervention2_GSA<- function(HBI_mean, betar, parous_mean, betam, betad,
                        sac_mean, xi, tau_mean,
                        N_b,coverage){
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
  #PA = exp(-(alpha_H_b*N_b + mu_A_b)*theta_d)###
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
  ## intervention
    chi = HBI_mean
    Pf = parous_mean
    A0 = sac_mean
    tau = tau_mean
    PA = 1 - A0
    betam = betam
    betad = betad
    pi = ifelse(betam+betad+betar<0.999, betam+betad+betar, 0.999)
    xi = xi
    #rho = rho
    c = coverage
    PA1 = A0*Pf*chi*eta*P_B_b*P_C_b/(P_B_b*P_C_b*(chi*P_D_b*P_E_b + (1-chi)*P_D_b*P_E_b)*eta*P_B_b*P_C_b)
    PAh = A0*Pf*(1-chi)/(eta*P_B_b*P_C_b*(chi*P_D_b*P_E_b +(1-chi)*P_D_b*P_E_b))
    PA2 = eta*PAh #encounter non-human host
    alpha_H_b = (1/N_b)*(PA1/(1-PA))*(-log(PA)/theta_d) #availability
    alpha_A_b = (1/N_b)*(PA2/(1-PA))*(-log(PA)/theta_d)
    mu_A_b = ((1-(PA + PA1 +PA2))/(1-PA))*(-log(PA)/theta_d)
    #PA+PA1+PA2+(1 - PA)*(mu_A_b/(((alpha_H_b+alpha_A_b)*N_b) + mu_A_b))###=1
    K_plus = floor(theta_s/tau) - 1
    #PA = exp(-((alpha_H_b+alpha_A_b)*N_b + mu_A_b)*theta_d)###
    #PAi = (1 - PA)*(((alpha_H_b+alpha_A_b)*N_b)/((alpha_H_b+alpha_A_b)*N_b + mu_A_b))
    #Pmu = (1 - PA)*(mu_A_b/((alpha_H_b+alpha_A_b)*N_b + mu_A_b))
    #Pdf = PAi*P_B_b*P_C_b*P_D_b*P_E_b
    MalariaHost = c(1,0,1,0,0)
    Hosts = which(MalariaHost == 1)
    N = c((1-c)*N_b, N_b, c*N_b, c*N_b,c*N_b)
    alpha_H = c(alpha_H_b, alpha_A_b, (1-pi)*alpha_H_b, betad*alpha_H_b, betam*alpha_H_b)
    P_C = c(P_C_b, P_C_b,(1-xi)*P_C_b, 1, 0)
    #mu_A = mu_A_b+rho*alpha_H_b*c*N_b
    PA = exp(-(sum(alpha_H*N) + mu_A_b)*theta_d)
    PAi = (1 - PA)*((alpha_H*N)/(sum(alpha_H*N) + mu_A_b))
    #PA + sum(PAi) + (1 - PA)*(mu_A_b/(sum(alpha_H*N) + mu_A_b)) ###=1 
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