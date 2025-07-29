################################################################################
#             SCRIPT PARA OS RESULTADOS INICIAIS DO PROCESSO INAR              #
################################################################################
rm(list = ls(all=T))
gc()
################################################################################
#                CARREGANDO AS BIBLIOTECAS NECESSÁRIAS                         #
################################################################################
library(optimx)
################################################################################
#                    PROCESSO POISSON INAR 1                                   #
################################################################################

simular_poisson <- function(n, alpha, lambda){
  x <- numeric(n)
  x[1] <- rpois(1, lambda) # Primeira obervação, ou valor inicial da série 
  
  for (t in 2:n){
    thinning <- rbinom(1, x[t-1], alpha) # Operador thinning binomial
    erro <- rpois(1, lambda) # parte do erro
    x[t] <- thinning + erro # Nova observação
  }
  return(x)
}

################################################################################
#                    ESTIMAÇÃO DE ALPHA E LAMBDA                               #
################################################################################

################################################################################
#           ESTIMAÇÃO ATRAVÉS DOS MINÍMOS QUADRADOS CONDICIONAIS               #
################################################################################


min_quadr_cond_poi <- function(x){
  n <- length(x)
  S1 <- sum(x[-1] * x[-n])  # Soma de X_t * X_{t-1}
  S2 <- sum(x[-n]^2)        # Soma de X_{t-1}^2
  S3 <- sum(x[-1])          # Soma de X_t
  S4 <- sum(x[-n])          # Soma de X_{t-1}
  
  # Estimador de aplha 
  
  alpha_hat <- (S1 - (S3 * S4) / (n-1)) / (S2 - (S4^2) / (n-1))
  alpha_hat <- max(0, min(alpha_hat, 1)) # Forçando para que alpha esteja entre 0 e 1 
  
  # Estimador de lambda 
  
  lambda_hat <- (S3 - alpha_hat * S4) / (n-1)
  lambda_hat <- max(0, lambda_hat) # Forçando para que lambda seja positivo
  
  return(c(alpha = alpha_hat, lamda = lambda_hat))
}

################################################################################
#     ESTIMAÇÃO ATRAVÉS DO ESTIMADOR DE MÁXIMA VEROSSIMILHANÇA CONDICIONAL     #
################################################################################

log_vero_poi <- function(parametros, x){
  alpha <- parametros[1]
  lambda <- parametros[2]
  n <- length(x)
  log_lik <- 0
  
  for(t in 2:n){
    probabilidade <- 0
    for(k in 0:min(x[t], x[t-1])){
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dpois(x[t] - k, lambda)
    }
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
} 

################################################################################
#             DEFININDO A ESTIMAÇÃO POR MÁXIMA VEROSSIMILHANÇA                 #
################################################################################

estimar_max_vero_poi <- function(x){
  # É preciso fornecer um chute inicial 
  
  chutes <- c(alpha = 0.5, lambda = mean(x))
  
  # Maximizando a função de log-verossimilhança 
  
  resultado <- optimx(
    par = chutes, 
    fn = log_vero_poi, 
    x = x, 
    method = "L-BFGS-B", 
    lower = c(0, 0), 
    upper = c(1, Inf)
  )
  
  return(resultado)
} 


################################################################################
#                             EXEMPLO DE USO                                   #
################################################################################

set.seed(2024)      # Setando a semente
n <- 200

alpha_verdadeiro <- 0.5
lambda_verdadeiro <- 4.0


################################################################################
#               PROCESSO INAR 1 GEOMÉTRICA PARAMETRIZADA PELA MÉDIA            #
################################################################################

simular_inar_geometrica <- function(n, alpha, mu){
  x <- numeric(n)
  p <- 1/ (1+ mu) # Fazendo a conversão de mu para o parâmetro p da geométrica 
  x[1] <- rgeom(1, p) # Primeira obeservação com distribuição geométrica 
  
  for(t in 2:n){
    thinning_binomial <- rbinom(1, x[t-1], alpha) # Operador Thinning binomial
    erro <- rgeom(1, p) # Termo do erro
    x[t] <- thinning_binomial + erro
  }
  return(x)
}

min_quadr_cond_geom <- function(x){
  n <- length(x)
  S1 <- sum(x[-1] * x[-n])  # Soma de X_t * X_{t-1}
  S2 <- sum(x[-n]^2)        # Soma de X_{t-1}^2
  S3 <- sum(x[-1])          # Soma de X_t
  S4 <- sum(x[-n])          # Soma de X_{t-1}
  
  # Estimador de aplha 
  
  alpha_hat <- (S1 - (S3 * S4) / (n-1)) / (S2 - (S4^2) / (n-1))
  alpha_hat <- max(0, min(alpha_hat, 1)) # Forçando para que alpha esteja entre 0 e 1 
  
  # Estimador de lambda 
  
  lambda_hat <- (S3 - alpha_hat * S4) / (n-1)
  lambda_hat <- max(0, lambda_hat) # Forçando para que lambda seja positivo
  
  return(c(alpha = alpha_hat, lamda = lambda_hat))
}

log_vero_geom <- function(parametros, x){
  alpha <- parametros[1]
  mu <- parametros[2]
  p <- 1/(1+mu)
  n <- length(x)
  log_lik <- 0
  
  for (t in 2:n) {
    probabilidade <- 0
    for (k in 0:min(x[t], x[t-1])) {
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dgeom(x[t] - k, p)
    }
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
}


estimar_max_vero_geom <- function(x){
  # Chutes iniciais 
  chutes <- c(alpha = 0.5, mu = mean(x))
  
  # maximizar a log-maxima verossimilhança 
  
  resultado <- optimx(
    par = chutes, 
    fn = log_vero_geom, 
    x = x, 
    method = "L-BFGS-B", 
    lower = c(0, 0), 
    upper = c(1, Inf)
  )
  return(resultado)
}



################################################################################
# PROCESSO INAR 1 COM DISTRIBUIÇÃO BINOMIAL NEGATIVA PARAMETRIZADA PELA MÉDIA  #
################################################################################

simular_inar_bin_negativa <- function(n, alpha, mu, sigma2){
  x <- numeric(n)
  r <- mu^2 / (sigma2 - mu)
  p <- mu/ sigma2
  x[1] <- rnbinom(1, size = r, prob = p) # Primeiro valor da simulação da série 
  
  for(t in 2:n){
    thinning_binomial <- rbinom(1, x[t-1], alpha) # Operador Thinning binomial
    erro <- rnbinom(1, size = r, prob = p) # Termo de erro
    x[t] <- thinning_binomial + erro # Valor atual
  }
  return(x)
}


min_quadr_cond_bin_negativa <- function(x){
  n <- length(x)
  S1 <- sum(x[-1] * x[-n])  # Soma de X_t * X_{t-1}
  S2 <- sum(x[-n]^2)        # Soma de X_{t-1}^2
  S3 <- sum(x[-1])          # Soma de X_t
  S4 <- sum(x[-n])          # Soma de X_{t-1}
  
  # Estimador de aplha 
  
  alpha_hat <- (S1 - (S3 * S4) / (n-1)) / (S2 - (S4^2) / (n-1))
  alpha_hat <- max(0, min(alpha_hat, 1)) # Forçando para que alpha esteja entre 0 e 1 
  
  # Estimador de lambda 
  
  lambda_hat <- (S3 - alpha_hat * S4) / (n-1)
  lambda_hat <- max(0, lambda_hat) # Forçando para que lambda seja positivo
  
  return(c(alpha = alpha_hat, lamda = lambda_hat))
}

log_vero_bin_negativa <- function(parametros, x){
  alpha <- parametros[1]
  mu <- parametros[2]
  sigma2 <- parametros[3]
  
  r <- mu^2 / (sigma2 - mu)
  p <- mu / sigma2
  n <- length(x)
  log_lik <- 0
  
  for (t in 2:n) {
    probabilidade <- 0
    for (k in 0:min(x[t], x[t-1])) {
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dnbinom(x[t] - k, size = r, prob = p)
    }
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
}

estimar_max_vero_bin_negativa <- function(x){
  # Chutes iniciais 
  chutes <- c(alpha = 0.5, mu = mean(x), sigma2 = var(x))
  
  resultado <- optimx(
    par = chutes, 
    fn = log_vero_bin_negativa, 
    x = x, 
    method = "L-GBFS-B", 
    lower = c(0, 0, 0), 
    upper = c(1, Inf, Inf)
  )
  return(resultado)
}
################################################################################
#                            SIMULANDO OS DADOS                                #
################################################################################

x <- simular_poisson(n, alpha = alpha_verdadeiro, lambda = lambda_verdadeiro)

resultados_minimos_quadrados_condicionais <- min_quadr_cond_poi(x)
print("Estimativas por Minímos Quadrados Condicionais: ")
print(resultados_minimos_quadrados_condicionais)

resultados_maxima_verossimilhança_condicionais <- estimar_max_vero_poi(x)

print("Estimativas por Máxima Verossimilhança Condicional: ")
print(resultados_maxima_verossimilhança_condicionais)

x_geom <- simular_inar_geometrica(n, alpha = alpha_verdadeiro, mu = lambda_verdadeiro)
resultados_minimos_quadrados_condicionais_geom <- min_quadr_cond_geom(x_geom)
paste0("Estimativas por Mínimos Quadrados Condicionais (Distribuição Geométrica): \n", 
       "Lambda: ", resultados_minimos_quadrados_condicionais_geom[1], " Mu: ", resultados_minimos_quadrados_condicionais_geom[2])


resultados_maxima_verossimilhança_condicionais_geom <- estimar_max_vero_geom(x_geom)
resultados_maxima_verossimilhança_condicionais_geom


x_binomial_negativa <- simular_inar_bin_negativa(n, alpha = alpha_verdadeiro, mu = 2, sigma2 = 4)
resultados_minimos_quadrados_condicionais_binomial_negativa <- min_quadr_cond_bin_negativa(x_binomial_negativa)
resultados_minimos_quadrados_condicionais_binomial_negativa

resultados_maxima_verossimilhança_condicionais_bin_negativa <- estimar_max_vero_bin_negativa(x_binomial_negativa)
resultados_maxima_verossimilhança_condicionais_bin_negativa
