################################################################################
#                    SCRIPT DAS SIMULAÇÕES PARA O TCC                          #
################################################################################
rm(list=ls(all=T))
gc()
################################################################################
#                        CARREGANDO AS BIBLIOTECAS                             #
################################################################################
library(optimx)
library(DataExplorer)
library(dygraphs)
library(plotly)
library(skimr)
library(gt)
library(ggplot2)
library(tidyverse)
library(xtable)
library(moments)
library(tsibble)
library(readxl)
library(openxlsx)
library(ggrepel)
library(forecast)
library(AER)
################################################################################
#                         CRIANDO AS FUNÇÕES                                   #
################################################################################
################################################################################
#                     SIMULAR CADA SÉRIE INDIVIDUAL                            #
################################################################################

simular_poisson <- function(n, alpha, lambda){
  x <- numeric(n)
  x[1] <- rpois(1, lambda)
  
  for(t in 2:n){
    thinnig_binomial <- rbinom(1, x[t-1], alpha)
    erro <- rpois(1, lambda)
    x[t] <- thinnig_binomial + erro
  }
  return(x)
}

simular_binomial_negativa <- function(n, alpha, mu, sigma2){
  x <- numeric(n)
  r <- mu^2 / (sigma2 - mu)
  p <- mu / sigma2
  x[1] <- rnbinom(1, size = r, prob = p)
  
  for (t in 2:n) {
    thinning_binomial <- rbinom(1, x[t-1], alpha)
    erro <- rnbinom(1, size = r, prob = p)
    x[t] <- thinning_binomial + erro
  }
  return(x)
}

simular_geometrica <- function(n, alpha, mu){
  x <- numeric(n)
  p <- 1 / (1 + mu)
  x[1] <- rgeom(1, p)
  
  for (t in 2:n) {
    thinning_binomial <- rbinom(1, x[t-1], alpha)
    erro <- rgeom(1, p)
    x[t] <- thinning_binomial + erro
  }
  return(x)
}
################################################################################
#                 FUNÇÕES DOS ESTIMADORES DE MQC E MVC                         #
################################################################################
minimos_quadrados <- function(x){
  n <- length(x)
  S1 <- sum(x[-1] * x[-n])                                  # Soma de Y_t * Y_{t-1}
  S2 <- sum(x[-n]^2)                                        # Soma de Y_{t-1}^2
  S3 <- sum(x[-1])                                          # Soma de Y_t
  S4 <- sum(x[-n])                                          # Soma de Y_{t-1}
  
  # Estimado de alpha 
  
  alpha_hat <- (S1 - (S3 * S4) / (n-1)) / (S2 - (S4^2) / (n-1))
  alpha_hat <- max(0, min(alpha_hat, 1)) # Por definição sabemos que alpha pertence ao intervalo [0,1]
  
  lambda_hat <- (S3 - alpha_hat * S4) / (n-1)
  lambda_hat <- max(0, lambda_hat) # Como valores de contagem são positivos, temos que definir esse intervalo
  
  return(c(alpha_estimado = alpha_hat, lambda_estimado = lambda_hat))
}
logveros_poisson <- function(parametros, x){
  if (any(is.na(parametros)) || any(is.infinite(parametros))) return(1e6)
  
  alpha <- parametros[1]
  lambda <- parametros[2]
  if (alpha < 0 || alpha > 1 || lambda <= 0) return(1e6)
  
  n <- length(x)
  log_lik <- 0
  
  for(t in 2:n){
    probabilidade <- 0
    for (k in 0:min(x[t], x[t-1])) {
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dpois(x[t] - k, lambda)
    }
    if (probabilidade <= 0 || is.na(probabilidade)) return(1e6)
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
}

estimar_max_ver_poisson <- function(x){
  chutes <- c(alpha = 0.4, lambda = mean(x))
  
  resultado <- optimx::optimx(
    par = chutes, 
    fn = logveros_poisson, 
    x = x, 
    method = "L-BFGS-B",
    lower = c(0, 0), 
    upper = c(1, Inf)
  )
  return(c(alpha = resultado$alpha, lambda = resultado$lambda))
}

logveros_binomial_negativa <- function(parametros, x){
  if (any(is.na(parametros)) || any(is.infinite(parametros))) return(1e6)
  
  alpha <- parametros[1] 
  mu <- parametros[2]
  sigma2 <- parametros[3]
  
  if (alpha < 0 || alpha > 1 || mu <= 0 || sigma2 <= mu) return(1e6)
  
  r <- mu^2 / (sigma2 - mu)
  p <- mu / sigma2
  n <- length(x)
  log_lik <- 0
  
  for (t in 2:n) {
    probabilidade <- 0
    for (k in 0:min(x[t], x[t-1])) {
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dnbinom(x[t] - k , size = r, prob = p)
    }
    if (probabilidade <= 0 || is.na(probabilidade)) return(1e6)
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
}

estimar_max_ver_binomial_negativa <- function(x){
  if (length(na.omit(x)) < 2) return(c(alpha = NA, mu = NA, sigma2 = NA))
  
  mu_inicial <- mean(x, na.rm = TRUE)
  sigma2_inicial <- var(x, na.rm = TRUE)
  
  if (is.na(mu_inicial) || is.na(sigma2_inicial) || sigma2_inicial <= mu_inicial) {
    return(c(alpha = NA, mu = NA, sigma2 = NA))
  }
  chutes <- c(alpha = 0.4, mu = mu_inicial, sigma2 = sigma2_inicial)
  
  resultado <- optimx::optimx(
    par = chutes, 
    fn = logveros_binomial_negativa, 
    x = x, 
    method = "L-BFGS-B", 
    lower = c(0, 0, 0), 
    upper = c(1, Inf, Inf)
  )
  return(c(alpha = resultado$alpha, mu = resultado$mu, sigma2 = resultado$sigma2))
}

logveros_geometrica <- function(parametros, x){
  if (any(is.na(parametros)) || any(is.infinite(parametros))) return(1e6)
  
  alpha <- parametros[1]
  mu <- parametros[2]
  
  if (alpha < 0 || alpha > 1 || mu <= 0) return(1e6)
  
  p <- 1 / (1 + mu)
  n <- length(x)
  log_lik <- 0
  
  for(t in 2:n){
    probabilidade <- 0
    for (k in 0:min(x[t], x[t-1])) {
      probabilidade <- probabilidade + dbinom(k, x[t-1], alpha) * dgeom(x[t] - k, p)
    }
    if (probabilidade <= 0 || is.na(probabilidade)) return(1e6)
    log_lik <- log_lik + log(probabilidade)
  }
  return(-log_lik)
}

estimar_max_ver_geometrica <- function(x){
  chutes <- c(alpha = 0.4, mu = mean(x))
  
  resultado <- optimx::optimx(
    par = chutes, 
    fn = logveros_geometrica, 
    x = x, 
    method = "L-BFGS-B", 
    lower = c(0, 0), 
    upper = c(1, Inf)
  )
  return(c(alpha = resultado$alpha, mu = resultado$mu))
}

################################################################################
#                     FUNÇÃO DE SIMULAÇÃO DE MONTE CARLO                       #
################################################################################

simulacao_monte_carlo_bin_negativa <- function(n_simulacoes, n_amostra, alpha_verdadeiro, mu_verdadeiro, sigma_verdadeiro){
  resultados_CLS <- matrix(NA, nrow = n_simulacoes, ncol = 2)
  resultados_CML <- matrix(NA, nrow = n_simulacoes, ncol = 3)
  
  
  for (i in 1:n_simulacoes) {
    x <- simular_binomial_negativa(n = n_amostra, alpha = alpha_verdadeiro, mu = mu_verdadeiro, sigma2 = sigma_verdadeiro)
    resultados_CLS[i, ] <- minimos_quadrados(x)
    resultados_CML[i, ] <- estimar_max_ver_binomial_negativa(x)
  }
  
  
  media_CLS <- colMeans(resultados_CLS, na.rm = T)
  variancia_CLS <- apply(resultados_CLS, 2, var, na.rm = T)
  vies_CLS <- media_CLS - c(alpha_verdadeiro, mu_verdadeiro)
  eqm_CLS <- variancia_CLS + vies_CLS^2
  
  media_CML <- colMeans(resultados_CML, na.rm = T)
  variancia_CML <- apply(resultados_CML, 2, var, na.rm = T)
  vies_CML <- media_CML - c(alpha_verdadeiro, mu_verdadeiro, sigma_verdadeiro)
  eqm_CML <- variancia_CML + vies_CML^2
  
  return(list(
    media_Minimos_Quadrados_condicionais = media_CLS, 
    variancia_Minimos_Quadrados_condicionais = variancia_CLS, 
    vies_minimos_quadrados_condicionais = vies_CLS,
    eqm_minimos_quadrados_condicionais = eqm_CLS,
    media_Maxima_verossimilhanca_condicional = media_CML, 
    variancia_Maxima_verossimilhanca_condicional = variancia_CML,
    vies_maxima_verossimilhanca_condicional = vies_CML, 
    eqm_maxima_verossimilhanca_condicional = eqm_CML,
    minimos_quadrados_resultado <- resultados_CLS, 
    maxima_verossimilhanca_resultado <- resultados_CML
  ))
}

simulacao_monte_carlo_poisson <- function(n_simulacoes, n_amostra, alpha_verdadeiro, lambda_verdadeiro){
  resultados_CLS <- matrix(NA, nrow = n_simulacoes, ncol = 2)
  resultados_CML <- matrix(NA, nrow = n_simulacoes, ncol = 2)
  
  
  for (i in 1:n_simulacoes) {
    x <- simular_poisson(n = n_amostra, alpha = alpha_verdadeiro, lambda = lambda_verdadeiro)
    resultados_CLS[i, ] <- minimos_quadrados(x)
    resultados_CML[i, ] <- estimar_max_ver_poisson(x)
  }
  
  
  media_CLS <- colMeans(resultados_CLS, na.rm = T)
  variancia_CLS <- apply(resultados_CLS, 2, var, na.rm = T)
  vies_CLS <- media_CLS - c(alpha_verdadeiro, lambda_verdadeiro)
  eqm_CLS <- variancia_CLS + vies_CLS^2
  
  media_CML <- colMeans(resultados_CML, na.rm = T)
  variancia_CML <- apply(resultados_CML, 2, var, na.rm = T)
  vies_CML <- media_CML - c(alpha_verdadeiro, lambda_verdadeiro)
  eqm_CML <- variancia_CML + vies_CML^2
  
  return(list(
    media_Minimos_Quadrados_condicionais = media_CLS, 
    variancia_Minimos_Quadrados_condicionais = variancia_CLS, 
    vies_minimos_quadrados_condicionais = vies_CLS,
    eqm_minimos_quadrados_condicionais = eqm_CLS,
    media_Maxima_verossimilhanca_condicional = media_CML, 
    variancia_Maxima_verossimilhanca_condicional = variancia_CML,
    vies_maxima_verossimilhanca_condicional = vies_CML, 
    eqm_maxima_verossimilhanca_condicional = eqm_CML,
    minimos_quadrados_resultado <- resultados_CLS, 
    maxima_verossimilhanca_resultado <- resultados_CML
  ))
}

simulacao_monte_carlo_geometrica <- function(n_simulacoes, n_amostra, alpha_verdadeiro, mu_verdadeiro, sigma_verdadeiro){
  resultados_CLS <- matrix(NA, nrow = n_simulacoes, ncol = 2)
  resultados_CML <- matrix(NA, nrow = n_simulacoes, ncol = 2)
  
  
  for (i in 1:n_simulacoes) {
    x <- simular_geometrica(n = n_amostra, alpha = alpha_verdadeiro, mu = mu_verdadeiro)
    resultados_CLS[i, ] <- minimos_quadrados(x)
    resultados_CML[i, ] <- estimar_max_ver_geometrica(x)
  }
  
  
  media_CLS <- colMeans(resultados_CLS, na.rm = T)
  variancia_CLS <- apply(resultados_CLS, 2, var, na.rm = T)
  vies_CLS <- media_CLS - c(alpha_verdadeiro, mu_verdadeiro)
  eqm_CLS <- variancia_CLS + vies_CLS^2
  
  media_CML <- colMeans(resultados_CML, na.rm = T)
  variancia_CML <- apply(resultados_CML, 2, var, na.rm = T)
  vies_CML <- media_CML - c(alpha_verdadeiro, mu_verdadeiro)
  eqm_CML <- variancia_CML + vies_CML^2
  
  return(list(
    media_Minimos_Quadrados_condicionais = media_CLS, 
    variancia_Minimos_Quadrados_condicionais = variancia_CLS, 
    vies_minimos_quadrados_condicionais = vies_CLS,
    eqm_minimos_quadrados_condicionais = eqm_CLS,
    media_Maxima_verossimilhanca_condicional = media_CML, 
    variancia_Maxima_verossimilhanca_condicional = variancia_CML,
    vies_maxima_verossimilhanca_condicional = vies_CML, 
    eqm_maxima_verossimilhanca_condicional = eqm_CML,
    minimos_quadrados_resultado <- resultados_CLS, 
    maxima_verossimilhanca_resultado <- resultados_CML
  ))
}

################################################################################
#                     FUNÇÕES DE PREVISÃO UM PASSO À FRENTE                   #
################################################################################

# Previsão para modelo Poisson
previsao_poisson <- function(x_anterior, alpha, lambda){
  return(alpha * x_anterior + lambda)
}

# Previsão para modelo Binomial Negativa
previsao_binomial_negativa <- function(x_anterior, alpha, mu, sigma2){
  return(alpha * x_anterior + mu)
}

# Previsão para modelo Geométrica
previsao_geometrica <- function(x_anterior, alpha, mu){
  return(alpha * x_anterior + mu)
}

# Função para calcular métricas de erro
calcular_metricas_erro <- function(valores_reais, valores_previstos){
  erro_absoluto <- abs(valores_reais - valores_previstos)
  erro_quadratico <- (valores_reais - valores_previstos)^2
  
  mae <- mean(erro_absoluto, na.rm = TRUE)
  rmse <- sqrt(mean(erro_quadratico, na.rm = TRUE))
  
  return(list(mae = mae, rmse = rmse))
}

################################################################################
#           SIMULAÇÃO MONTE CARLO PARA PREVISÃO UM PASSO À FRENTE             #
################################################################################

# Simulação Monte Carlo para previsão - Modelo Poisson
simulacao_monte_carlo_previsao_poisson <- function(n_simulacoes, n_treino, n_teste, alpha_verdadeiro, lambda_verdadeiro){
  mae_CLS <- numeric(n_simulacoes)
  rmse_CLS <- numeric(n_simulacoes)
  mae_CML <- numeric(n_simulacoes)
  rmse_CML <- numeric(n_simulacoes)
  
  for (i in 1:n_simulacoes) {
    # Gerar série temporal completa
    x_completa <- simular_poisson(n = n_treino + n_teste, alpha = alpha_verdadeiro, lambda = lambda_verdadeiro)
    
    # Dividir em treino e teste
    x_treino <- x_completa[1:n_treino]
    x_teste <- x_completa[(n_treino + 1):(n_treino + n_teste)]
    
    # Estimar parâmetros nos dados de treino
    parametros_CLS <- minimos_quadrados(x_treino)
    parametros_CML <- estimar_max_ver_poisson(x_treino)
    
    # Fazer previsões um passo à frente
    previsoes_CLS <- numeric(n_teste)
    previsoes_CML <- numeric(n_teste)
    
    for (t in 1:n_teste) {
      x_anterior <- ifelse(t == 1, x_treino[n_treino], x_completa[n_treino + t - 1])
      
      previsoes_CLS[t] <- previsao_poisson(x_anterior, parametros_CLS[1], parametros_CLS[2])
      previsoes_CML[t] <- previsao_poisson(x_anterior, parametros_CML[1], parametros_CML[2])
    }
    
    # Calcular métricas de erro
    metricas_CLS <- calcular_metricas_erro(x_teste, previsoes_CLS)
    metricas_CML <- calcular_metricas_erro(x_teste, previsoes_CML)
    
    mae_CLS[i] <- metricas_CLS$mae
    rmse_CLS[i] <- metricas_CLS$rmse
    mae_CML[i] <- metricas_CML$mae
    rmse_CML[i] <- metricas_CML$rmse
  }
  
  return(list(
    mae_media_CLS = mean(mae_CLS, na.rm = TRUE),
    rmse_media_CLS = mean(rmse_CLS, na.rm = TRUE),
    mae_media_CML = mean(mae_CML, na.rm = TRUE),
    rmse_media_CML = mean(rmse_CML, na.rm = TRUE),
    mae_todos_CLS = mae_CLS,
    rmse_todos_CLS = rmse_CLS,
    mae_todos_CML = mae_CML,
    rmse_todos_CML = rmse_CML
  ))
}

# Simulação Monte Carlo para previsão - Modelo Binomial Negativa
simulacao_monte_carlo_previsao_bin_negativa <- function(n_simulacoes, n_treino, n_teste, alpha_verdadeiro, mu_verdadeiro, sigma_verdadeiro){
  mae_CLS <- numeric(n_simulacoes)
  rmse_CLS <- numeric(n_simulacoes)
  mae_CML <- numeric(n_simulacoes)
  rmse_CML <- numeric(n_simulacoes)
  
  for (i in 1:n_simulacoes) {
    # Gerar série temporal completa
    x_completa <- simular_binomial_negativa(n = n_treino + n_teste, alpha = alpha_verdadeiro, mu = mu_verdadeiro, sigma2 = sigma_verdadeiro)
    
    # Dividir em treino e teste
    x_treino <- x_completa[1:n_treino]
    x_teste <- x_completa[(n_treino + 1):(n_treino + n_teste)]
    
    # Estimar parâmetros nos dados de treino
    parametros_CLS <- minimos_quadrados(x_treino)
    parametros_CML <- estimar_max_ver_binomial_negativa(x_treino)
    
    # Fazer previsões um passo à frente
    previsoes_CLS <- numeric(n_teste)
    previsoes_CML <- numeric(n_teste)
    
    for (t in 1:n_teste) {
      x_anterior <- ifelse(t == 1, x_treino[n_treino], x_completa[n_treino + t - 1])
      
      previsoes_CLS[t] <- previsao_binomial_negativa(x_anterior, parametros_CLS[1], mu_verdadeiro)
      if (!any(is.na(parametros_CML))) {
        previsoes_CML[t] <- previsao_binomial_negativa(x_anterior, parametros_CML[1], parametros_CML[2])
      } else {
        previsoes_CML[t] <- NA
      }
    }
    
    # Calcular métricas de erro
    metricas_CLS <- calcular_metricas_erro(x_teste, previsoes_CLS)
    metricas_CML <- calcular_metricas_erro(x_teste, previsoes_CML)
    
    mae_CLS[i] <- metricas_CLS$mae
    rmse_CLS[i] <- metricas_CLS$rmse
    mae_CML[i] <- metricas_CML$mae
    rmse_CML[i] <- metricas_CML$rmse
  }
  
  return(list(
    mae_media_CLS = mean(mae_CLS, na.rm = TRUE),
    rmse_media_CLS = mean(rmse_CLS, na.rm = TRUE),
    mae_media_CML = mean(mae_CML, na.rm = TRUE),
    rmse_media_CML = mean(rmse_CML, na.rm = TRUE),
    mae_todos_CLS = mae_CLS,
    rmse_todos_CLS = rmse_CLS,
    mae_todos_CML = mae_CML,
    rmse_todos_CML = rmse_CML
  ))
}

# Simulação Monte Carlo para previsão - Modelo Geométrica
simulacao_monte_carlo_previsao_geometrica <- function(n_simulacoes, n_treino, n_teste, alpha_verdadeiro, mu_verdadeiro){
  mae_CLS <- numeric(n_simulacoes)
  rmse_CLS <- numeric(n_simulacoes)
  mae_CML <- numeric(n_simulacoes)
  rmse_CML <- numeric(n_simulacoes)
  
  for (i in 1:n_simulacoes) {
    # Gerar série temporal completa
    x_completa <- simular_geometrica(n = n_treino + n_teste, alpha = alpha_verdadeiro, mu = mu_verdadeiro)
    
    # Dividir em treino e teste
    x_treino <- x_completa[1:n_treino]
    x_teste <- x_completa[(n_treino + 1):(n_treino + n_teste)]
    
    # Estimar parâmetros nos dados de treino
    parametros_CLS <- minimos_quadrados(x_treino)
    parametros_CML <- estimar_max_ver_geometrica(x_treino)
    
    # Fazer previsões um passo à frente
    previsoes_CLS <- numeric(n_teste)
    previsoes_CML <- numeric(n_teste)
    
    for (t in 1:n_teste) {
      x_anterior <- ifelse(t == 1, x_treino[n_treino], x_completa[n_treino + t - 1])
      
      previsoes_CLS[t] <- previsao_geometrica(x_anterior, parametros_CLS[1], mu_verdadeiro)
      previsoes_CML[t] <- previsao_geometrica(x_anterior, parametros_CML[1], parametros_CML[2])
    }
    
    # Calcular métricas de erro
    metricas_CLS <- calcular_metricas_erro(x_teste, previsoes_CLS)
    metricas_CML <- calcular_metricas_erro(x_teste, previsoes_CML)
    
    mae_CLS[i] <- metricas_CLS$mae
    rmse_CLS[i] <- metricas_CLS$rmse
    mae_CML[i] <- metricas_CML$mae
    rmse_CML[i] <- metricas_CML$rmse
  }
  
  return(list(
    mae_media_CLS = mean(mae_CLS, na.rm = TRUE),
    rmse_media_CLS = mean(rmse_CLS, na.rm = TRUE),
    mae_media_CML = mean(mae_CML, na.rm = TRUE),
    rmse_media_CML = mean(rmse_CML, na.rm = TRUE),
    mae_todos_CLS = mae_CLS,
    rmse_todos_CLS = rmse_CLS,
    mae_todos_CML = mae_CML,
    rmse_todos_CML = rmse_CML
  ))
}
################################################################################
#                               SIMULAÇOES                                     #
################################################################################
set.seed(1234) # Definindo a semente 

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - POISSON - N = 50                 #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
  Estimador = rep(c("alpha", "mu"), times = 2),
  Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
  Média = c(resultados$media_Minimos_Quadrados_condicionais,
            resultados$media_Maxima_verossimilhanca_condicional),
  Viés = c(resultados$vies_minimos_quadrados_condicionais,
           resultados$vies_maxima_verossimilhanca_condicional),
  Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                resultados$variancia_Maxima_verossimilhanca_condicional),
  EQM = c(resultados$eqm_minimos_quadrados_condicionais,
          resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.3 e lambda 1.5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - POISSON - N = 100                #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.3 e lambda 1.5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - POISSON - N = 300                #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.3 e lambda 1.5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - POISSON - N = 500                #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.3 e lambda 1.5, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - GEOMETRICA - N = 50              #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.3 e lambda 1.5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - GEOMETRICA - N = 100             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.3 e lambda 1.5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - GEOMETRICA - N = 300             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.3 e lambda 1.5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - GEOMETRICA - N = 500               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.3 e lambda 1.5, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3, LAMBDA 1.5 E SIGMA 2 - BIN NEGATIVA - N = 50   #
################################################################################
alpha_verdadeiro <- 0.3; lambda_verdadeiro <- 1.5; sigma_verdadeiro <- 2
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.3, lambda 1.5 e sigma 2, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - BIN NEGATIVA - N = 100           #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.3, lambda 1.5 e sigma 2, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - BIN NEGATIVA - N = 300           #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.3, lambda 1.5 e sigma 2, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.3 E LAMBDA 1.5 - BIN NEGATIVA - N = 500           #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.3, lambda 1.5 e sigma 2, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - POISSON - N = 50                   #
################################################################################
alpha_verdadeiro <- 0.5; lambda_verdadeiro <- 2
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.5 e lambda 2, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - POISSON - N = 100                  #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.5 e lambda 2, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - POISSON - N = 300                  #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.5 e lambda 2, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - POISSON - N = 500                  #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.3 e lambda 8, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - GEOMETRICA - N = 50                #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.5 e lambda 2, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - GEOMETRICA - N = 100               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.5 e lambda 2, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.3 E LAMBDA 8 - GEOMETRICA - N = 300               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.5 e lambda 2, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5 E LAMBDA 2 - GEOMETRICA - N = 500               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.5 e lambda 2, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5, LAMBDA 2 E SIGMA 3 - BIN NEGATIVA - N = 50     #
################################################################################
alpha_verdadeiro <- 0.5; lambda_verdadeiro <- 2; sigma_verdadeiro <- 3
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.5, lambda 2 e sigma 3, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.5, LAMBDA 2 E SIGMA 3 - BIN NEGATIVA - N = 100    #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.5, lambda 2 e sigma 3, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.5, LAMBDA 2 E SIGMA 3 - BIN NEGATIVA - N = 300    #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.5, lambda 2 e sigma 3, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.5, LAMBDA 2 E SIGMA 3 - BIN NEGATIVA - N = 500    #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.5, lambda 2 e sigma 3, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                  ALPHA 0.7 E LAMBDA 2.5 - POISSON - N = 50                   #
################################################################################
alpha_verdadeiro <- 0.7; lambda_verdadeiro <- 2.5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.7 e lambda 2.5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.7 E LAMBDA 2.5 - POISSON - N = 100                #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.7 e lambda 2.5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.7 E LAMBDA 2.5 - POISSON - N = 300                #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.7 e lambda 2.5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                  ALPHA 0.7 E LAMBDA 2.5 - POISSON - N = 500                  #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.7 e lambda 2.5, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.7 E LAMBDA 2.5 - GEOMETRICA - N = 50              #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.7 e lambda 2.5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                  ALPHA 0.7 E LAMBDA 2.5 - GEOMETRICA - N = 100               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.7 e lambda 2.5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.7 E LAMBDA 2.5 - GEOMETRICA - N = 300             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.7 e lambda 2.5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.7 E LAMBDA 2.5 - GEOMETRICA - N = 500             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.7 e lambda 2.5, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.7, LAMBDA 2.5 E SIGMA 4 - BIN NEGATIVA - N = 50   #
################################################################################
alpha_verdadeiro <- 0.7; lambda_verdadeiro <- 2.5; sigma_verdadeiro <- 4
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.7, lambda 2.5 e sigma 4, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                  ALPHA 0.7, LAMBDA 2.5 E SIGMA 4  - BIN NEGATIVA - N = 100   #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.7, lambda 2.5 e sigma 4, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                   ALPHA 0.7, LAMBDA 2.5 E SIGMA 4  - BIN NEGATIVA - N = 300  #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.7, lambda 2.5 e sigma 4, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                  ALPHA 0.7, LAMBDA 2.5 E SIGMA 4   - BIN NEGATIVA - N = 500  #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.7, lambda 2.5 e sigma 4, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

Sys.time()
################################################################################
#                    ALPHA 0.1 E LAMBDA 5 - POISSON - N = 50                 #
################################################################################
tempo_inicial <- paste0(format(Sys.time(), "%d/%m/%Y"), " em ", format(Sys.time(), "%H:%M:%S"))

alpha_verdadeiro <- 0.1; lambda_verdadeiro <- 5
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.1 e lambda 5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5 - - POISSON - N = 100              #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.1 e lambda 5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5   - - POISSON - N = 300              #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.1 e lambda 5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5   - POISSON - N = 500                #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.1 e lambda 5, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5 - GEOMETRICA - N = 50              #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.1 e lambda 5, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5   - GEOMETRICA - N = 100             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.1 e lambda 5, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.1 E LAMBDA 5 - GEOMETRICA - N = 300             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.1 e lambda 5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1 E LAMBDA 5 - GEOMETRICA - N = 500             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.7 e lambda 0.3, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.1, LAMBDA 5 E SIGMA 7     - BIN NEGATIVA - N = 50 #
################################################################################
alpha_verdadeiro <- 0.1; lambda_verdadeiro <- 5; sigma_verdadeiro <- 7
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.1, lambda 5 e sigma 7, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                ALPHA 0.1, LAMBDA 5 E SIGMA 7      - BIN NEGATIVA - N = 100   #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.1, lambda 5 e sigma 7, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                 ALPHA 0.1, LAMBDA 5 E SIGMA 7  - BIN NEGATIVA - N = 300      #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.7, lambda 0.3 e sigma 0.5, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                ALPHA 0.1, LAMBDA 5 E SIGMA 7       - BIN NEGATIVA - N = 500  #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.1, lambda 5 e sigma 7, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

tempo_final <- paste0(format(Sys.time(), "%d/%m/%Y"), " em ", format(Sys.time(), "%H:%M:%S"))

################################################################################
#                    ALPHA 0.9 E LAMBDA 4 - POISSON - N = 50                 #
################################################################################
tempo_inicial <- paste0(format(Sys.time(), "%d/%m/%Y"), " em ", format(Sys.time(), "%H:%M:%S"))

alpha_verdadeiro <- 0.9; lambda_verdadeiro <- 4
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.9 e lambda 4, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4 - - POISSON - N = 100                #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.9 e lambda 4, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4   - - POISSON - N = 300              #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.9 e lambda 4, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4   - POISSON - N = 500                #
################################################################################
resultados <- simulacao_monte_carlo_poisson(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, lambda_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Poisson - alpha 0.9 e lambda 4, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4 - GEOMETRICA - N = 50                #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.9 e lambda 4, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4   - GEOMETRICA - N = 100             #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.9 e lambda 4, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                    ALPHA 0.9 E LAMBDA 4 - GEOMETRICA - N = 300               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.9 e lambda 4, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9 E LAMBDA 4 - GEOMETRICA - N = 500               #
################################################################################
resultados <- simulacao_monte_carlo_geometrica(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Geometrica - alpha 0.9 e lambda 4, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                    ALPHA 0.9, LAMBDA 4 E SIGMA 4.1     - BIN NEGATIVA - N = 50 #
################################################################################
alpha_verdadeiro <- 0.9; lambda_verdadeiro <- 4; sigma_verdadeiro <- 4.1
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 50, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.9, lambda 4 e sigma 4.1, 50 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

################################################################################
#                ALPHA 0.9, LAMBDA 4 E SIGMA 4.1      - BIN NEGATIVA - N = 100   #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 100, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.9, lambda 4 e sigma 4.1, 100 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                 ALPHA 0.9, LAMBDA 4 E SIGMA 4.1 - BIN NEGATIVA - N = 300     #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 300, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.9, lambda 4 e sigma 4.1, 300 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt
################################################################################
#                ALPHA 0.9, LAMBDA 4 E SIGMA 4.1     - BIN NEGATIVA - N = 500  #
################################################################################
resultados <- simulacao_monte_carlo_bin_negativa(n_simulacoes = 5000, n_amostra = 500, alpha_verdadeiro = alpha_verdadeiro, mu_verdadeiro = lambda_verdadeiro, sigma_verdadeiro = sigma_verdadeiro)

tabela_resultados <- 
  data.frame(
    Estimador = rep(c("alpha", "mu"), times = 2),
    Método = rep(c("MQC (CLS)", "MVC (CML)"), each = 2),
    Média = c(resultados$media_Minimos_Quadrados_condicionais,
              resultados$media_Maxima_verossimilhanca_condicional[1:2]),
    Viés = c(resultados$vies_minimos_quadrados_condicionais,
             resultados$vies_maxima_verossimilhanca_condicional[1:2]),
    Variância = c(resultados$variancia_Minimos_Quadrados_condicionais,
                  resultados$variancia_Maxima_verossimilhanca_condicional[1:2]),
    EQM = c(resultados$eqm_minimos_quadrados_condicionais,
            resultados$eqm_maxima_verossimilhanca_condicional[1:2]))



tabela_resultados <- tabela_resultados[!(tabela_resultados$Estimador == "sigma²" & tabela_resultados$Método == "MQC (CLS)"), ]

tabela_gt <- tabela_resultados |>
  gt() |>
  tab_header(
    title = "Resumo dos Resultados das Simulações de Monte Carlo - Binomial Negativa - alpha 0.9, lambda 4 e sigma 4.1, 500 Observações",
    subtitle = "Média, Viés, Variância e EQM dos Estimadores"
  ) |>
  fmt_number(columns = c(Média, Viés, Variância, EQM), decimals = 4) |>
  cols_label(
    Estimador = "Parâmetro",
    Método = "Método",
    Média = "Média",
    Viés = "Viés",
    Variância = "Variância",
    EQM = "EQM"
  )

tabela_gt

tempo_final <- paste0(format(Sys.time(), "%d/%m/%Y"), " em ", format(Sys.time(), "%H:%M:%S"))

################################################################################
#                                 APLICAÇÃO PRÁTICA                            #
################################################################################

banco <- read_excel(file.choose())

#banco <- data.frame(data = seq.Date(from = as.Date("01-01-2025", "%d-%m-%Y"), to = as.Date("31-12-2025", "%d-%m-%Y"), by = "day"), quantidade = sample(1:34, 365, replace = T))

dados <- banco %>%
  mutate(data = as.Date(Data)) %>%
  mutate(quantidade = Quantidade) %>% 
  select(-c(Data, Quantidade)) %>% 
  arrange(data)



serie_temporal <- ts(dados$quantidade, 
                     start = c(year(min(dados$data)), month(min(dados$data)), day(min(dados$data))), 
                     frequency = 12)

serie_temporal
################################################################################
#                            ESTATÍSTICA DESCRITIVA                            #
################################################################################

estatisticas <- data.frame(
  Estatística = c("Obersavções", "Média", "Mediana", "Desvio Padrão", "Coeficiente de Variação", "Mínimo", "Máximo", "Assimetria", "Curtose"), 
  Valor = c(
    as.integer(length(serie_temporal)), 
    round(mean(serie_temporal, na.rm = T),3), 
    round(median(serie_temporal, na.rm = T),3), 
    round(sd(serie_temporal, na.rm = T),3), 
    round((sd(serie_temporal, na.rm = T)/mean(serie_temporal, na.rm = T)) * 100, 4), 
    round(min(serie_temporal, na.rm = T),3), 
    round(max(serie_temporal, na.rm = T),3), 
    round(skewness(serie_temporal, na.rm = T),3), 
    round(kurtosis(serie_temporal, na.rm = T),3)
  )
)

print(estatisticas)

tabela_latex <- xtable(estatisticas, 
                       caption = "Estatísticas Descritivas da Série de Seguros Cancelados por Dia", 
                       label = "tab:estatisticas_descritivas_seguros_cancelados")

print(tabela_latex, 
      include.rownames = FALSE)

ponto_maximo <- dados %>%
  filter(quantidade == max(quantidade)) %>%
  mutate(rotulo = paste0("Máximo: ", quantidade, "\n", format(data, "%d/%m/%Y")))

grafico_serie <- ggplot(data = dados, aes(x = data, y = quantidade)) +
  geom_line(color = "darkblue", linewidth = 0.7, alpha = 0.9) +
  geom_point(data = ponto_maximo, color = "#e74c3c", size = 2.5) +
  geom_label_repel(
    data = ponto_maximo,
    aes(label = rotulo),
    color = "#e74c3c",
    size = 3,
    nudge_y = max(dados$quantidade, na.rm = TRUE) * 0.08,
    segment.color = "#e74c3c",
    segment.size = 0.5,
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  geom_hline(yintercept = mean(dados$quantidade, na.rm = TRUE), 
             color = "#2c3e50", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = min(dados$data), 
           y = mean(dados$quantidade, na.rm = TRUE) * 1.02, 
           label = paste0("Média: ", round(mean(dados$quantidade, na.rm = TRUE), 2)), 
           hjust = 0, vjust = 0, color = "darkred", size = 3) +
  labs(
    title = "",
    subtitle = "", 
    x = "Período",
    y = "Quantidade",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(grafico_serie)

ggsave("serie_temporal_diaria_inar.png", 
       plot = grafico_serie, 
       width = 12,
       height = 6, 
       dpi = 300)

################################################################################
#                             FAZENDO A MODELAGEM                              #
################################################################################

cat("\n1. ESTIMAÇÃO POR MÍNIMOS QUADRADOS CONDICIONAIS (MQC):\n")
resultados_mqc <- minimos_quadrados(serie_temporal)
cat("Alpha estimado (MQC):", round(resultados_mqc[1], 4), "\n")
cat("Lambda estimado (MQC):", round(resultados_mqc[2], 4), "\n")

cat("\n2. ESTIMAÇÃO POR MÁXIMA VEROSSIMILHANÇA CONDICIONAL (MVC - Poisson):\n")
resultados_mvc_poisson <- estimar_max_ver_poisson(serie_temporal)
cat("Alpha estimado (MVC Poisson):", round(resultados_mvc_poisson[1], 4), "\n")
cat("Lambda estimado (MVC Poisson):", round(resultados_mvc_poisson[2], 4), "\n")


cat("\n3. ESTIMAÇÃO POR MÁXIMA VEROSSIMILHANÇA CONDICIONAL (MVC - Binomial Negativa):\n")
resultados_mvc_binneg <- estimar_max_ver_binomial_negativa(serie_temporal)
cat("Alpha estimado (MVC Bin Neg):", round(resultados_mvc_binneg[1], 4), "\n")
cat("Mu estimado (MVC Bin Neg):", round(resultados_mvc_binneg[2], 4), "\n")
cat("Sigma2 estimado (MVC Bin Neg):", round(resultados_mvc_binneg[3], 4), "\n")


cat("\n4. ESTIMAÇÃO POR MÁXIMA VEROSSIMILHANÇA CONDICIONAL (MVC - Geométrica):\n")
resultados_mvc_geo <- estimar_max_ver_geometrica(serie_temporal)
cat("Alpha estimado (MVC Geo):", round(resultados_mvc_geo[1], 4), "\n")
cat("Mu estimado (MVC Geo):", round(resultados_mvc_geo[2], 4), "\n")


calcular_aic <- function(log_ver, k) {
  return(2 * k - 2 * log_ver)
}

calcular_bic <- function(log_ver, k, n) {
  return(log(n) * k - 2 * log_ver)
}

n <- length(serie_temporal)
log_ver_poisson <- -logveros_poisson(resultados_mvc_poisson, serie_temporal)
log_ver_binneg <- -logveros_binomial_negativa(resultados_mvc_binneg, serie_temporal)
log_ver_geo <- -logveros_geometrica(resultados_mvc_geo, serie_temporal)

aic_poisson <- calcular_aic(log_ver_poisson, 2)
aic_binneg <- calcular_aic(log_ver_binneg, 3)
aic_geo <- calcular_aic(log_ver_geo, 2)

bic_poisson <- calcular_bic(log_ver_poisson, 2, n)
bic_binneg <- calcular_bic(log_ver_binneg, 3, n)
bic_geo <- calcular_bic(log_ver_geo, 2, n)

cat("\n=== CRITÉRIOS DE INFORMAÇÃO PARA SELEÇÃO DE MODELOS ===\n")
cat("Poisson - AIC:", round(aic_poisson, 2), "BIC:", round(bic_poisson, 2), "\n")
cat("Binomial Negativa - AIC:", round(aic_binneg, 2), "BIC:", round(bic_binneg, 2), "\n")
cat("Geométrica - AIC:", round(aic_geo, 2), "BIC:", round(bic_geo, 2), "\n")


modelos <- c("Poisson", "Binomial Negativa", "Geométrica")
aics <- c(aic_poisson, aic_binneg, aic_geo)
bics <- c(bic_poisson, bic_binneg, bic_geo)

melhor_aic <- modelos[which.min(aics)]
melhor_bic <- modelos[which.min(bics)]
melhor_modelo <- modelos[which.min(aics)]

cat("Melhor modelo:", melhor_modelo, "\n")
cat("Melhor modelo por AIC:", melhor_aic, "\n")
cat("Melhor modelo por BIC:", melhor_bic, "\n")

if (melhor_modelo == "Poisson") {
  alpha <- resultados_mvc_poisson[1]
  lambda <- resultados_mvc_poisson[2]
  mu_hat <- sapply(2:length(serie_temporal), function(t) round(alpha * serie_temporal[t-1] + lambda), 0)
  
} else if (melhor_modelo == "Binomial Negativa") {
  alpha <- resultados_mvc_binneg[1]
  mu <- resultados_mvc_binneg[2]
  mu_hat <- sapply(2:length(serie_temporal), function(t) round(alpha * serie_temporal[t-1] + mu), 0)
  
} else {
  alpha <- resultados_mvc_geo[1]
  mu <- resultados_mvc_geo[2]
  mu_hat <- sapply(2:length(serie_temporal), function(t) round(alpha * serie_temporal[t-1] + mu, 0))
}


y <- serie_temporal[2:length(serie_temporal)]
residuos_pearson <- (y - mu_hat) / sqrt(mu_hat)


df_res <- data.frame(
  tempo = 1:length(residuos_pearson),
  residuos = residuos_pearson
)

p1 <- ggplot(df_res, aes(x = tempo, y = residuos)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "", x = "Observação", y = "Resíduo") +
  theme_classic()
print(p1)

ggsave("residuos_pearson.png", p1, width = 7, height = 4, dpi = 300)


p2 <- ggAcf(residuos_pearson) + 
  labs(title = "") +
  theme_classic(base_size = 14)

print(p2)

ggsave("acf_residuos_pearson.png", p2, width = 7, height = 4, dpi = 300)


ljung <- Box.test(residuos_pearson, lag = 10, type = "Ljung-Box")
cat("\n=== TESTE DE RESÍDUOS ===\n")
print(ljung)

Box.test(serie_temporal, lag = 10, type = "Ljung-Box")

if (melhor_aic == "Poisson") {
  cat("\n=== PREVISÕES COM MODELO POISSON ===\n")
  alpha <- resultados_mvc_poisson[1]
  lambda <- resultados_mvc_poisson[2]
  previsoes <- previsao_poisson(tail(serie_temporal, 1), alpha, lambda)
  cat("Previsão 1 passo à frente:", round(previsoes, 0), "\n")
  
} else if (melhor_aic == "Binomial Negativa") {
  cat("\n=== PREVISÕES COM MODELO BINOMIAL NEGATIVA ===\n")
  alpha <- resultados_mvc_binneg[1]
  mu <- resultados_mvc_binneg[2]
  previsoes <- previsao_binomial_negativa(tail(serie_temporal, 1), alpha, mu)
  cat("Previsão 1 passo à frente:", round(previsoes, 0), "\n")
  
} else {
  cat("\n=== PREVISÕES COM MODELO GEOMÉTRICA ===\n")
  alpha <- resultados_mvc_geo[1]
  mu <- resultados_mvc_geo[2]
  previsoes <- previsao_geometrica(tail(serie_temporal, 1), alpha, mu)
  cat("Previsão 1 passo à frente:", round(previsoes, 0), "\n")
}

horizontes <- c(1, 3, 6, 12)

reestimar_parametros <- function(serie, modelo, h) {
  serie_treino <- serie[1:(length(serie) - h)]
  
  parametros <- NULL
  
  if (modelo == "Poisson") {
    if (exists("estimar_max_ver_poisson")) {
      parametros <- estimar_max_ver_poisson(serie_treino)
    } else {
      stop("Função estimar_max_ver_poisson não encontrada")
    }
  } else if (modelo == "Binomial_Negativa") {
    if (exists("estimar_max_ver_binomial_negativa")) {
      parametros <- estimar_max_ver_binomial_negativa(serie_treino)
    } else {
      stop("Função estimar_max_ver_binomial_negativa não encontrada")
    }
  } else if (modelo == "Geométrica") {
    if (exists("estimar_max_ver_geometrica")) {
      parametros <- estimar_max_ver_geometrica(serie_treino)
    } else {
      stop("Função estimar_max_ver_geometrica não encontrada")
    }
  } else {
    stop("Modelo não reconhecido: ", modelo)
  }
  
  if (is.null(parametros)) {
    stop("Falha na estimação dos parâmetros para o modelo: ", modelo)
  }
  
  return(parametros)
}


previsao_poisson_multi <- function(x_inicial, alpha, lambda, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + lambda, 0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}

previsao_binomial_negativa_multi <- function(x_inicial, alpha, mu, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + mu, 0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}

previsao_geometrica_multi <- function(x_inicial, alpha, mu, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + mu, 0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}


previsoes_multi <- list()
parametros_reestimados <- list()
valores_reais <- list()
residuos <- list()

cat("\n=== VALIDAÇÃO WALK-FORWARD PARA TODOS OS HORIZONTES ===\n")

for (h in horizontes) {
  cat("\n--- Horizonte:", h, "passo(s) ---\n")
  
  if (h <= length(serie_temporal)) {

    parametros_reestimados[[paste0("h", h)]] <- reestimar_parametros(serie_temporal, melhor_modelo, h)
    cat("Parâmetros reestimados:", round(parametros_reestimados[[paste0("h", h)]], 4), "\n")
    
    valor_inicial <- tail(serie_temporal[1:(length(serie_temporal) - h)], 1)
    cat("Valor inicial para previsão:", valor_inicial, "\n")
    
    if (melhor_modelo == "Poisson") {
      previsoes_multi[[paste0("h", h)]] <- previsao_poisson_multi(valor_inicial, 
                                                                  parametros_reestimados[[paste0("h", h)]][1], 
                                                                  parametros_reestimados[[paste0("h", h)]][2], h)
    } else if (melhor_modelo == "Binomial_Negativa") {
      previsoes_multi[[paste0("h", h)]] <- previsao_binomial_negativa_multi(valor_inicial, 
                                                                            parametros_reestimados[[paste0("h", h)]][1], 
                                                                            parametros_reestimados[[paste0("h", h)]][2], h)
    } else {
      previsoes_multi[[paste0("h", h)]] <- previsao_geometrica_multi(valor_inicial, 
                                                                     parametros_reestimados[[paste0("h", h)]][1], 
                                                                     parametros_reestimados[[paste0("h", h)]][2], h)
    }
    
    
    valores_reais[[paste0("h", h)]] <- tail(serie_temporal, h)
    
    
    residuos[[paste0("h", h)]] <- valores_reais[[paste0("h", h)]] - previsoes_multi[[paste0("h", h)]]
    
    cat("Previsão:", previsoes_multi[[paste0("h", h)]], "\n")
    cat("Valor real:", valores_reais[[paste0("h", h)]], "\n")
    cat("Resíduos:", residuos[[paste0("h", h)]], "\n")
  }
}


cat("\n=== ANÁLISE DETALHADA POR HORIZONTE ===\n")

datas_reais <- dados$data  

for (h in horizontes) {
  if (h <= length(serie_temporal)) {
    cat("\n--- Análise para horizonte", h, "---\n")
    
    mae <- mean(abs(residuos[[paste0("h", h)]]), na.rm = TRUE)
    rmse <- sqrt(mean(residuos[[paste0("h", h)]]^2, na.rm = TRUE))
    mape <- mean(abs(residuos[[paste0("h", h)]]/valores_reais[[paste0("h", h)]]), na.rm = TRUE) * 100
    
    cat("MAE:", round(mae, 4), "\n")
    cat("RMSE:", round(rmse, 4), "\n")
    cat("MAPE:", round(mape, 2), "%\n")
    
    
    datas_previsao <- tail(datas_reais, h)  
    
  
    dados_comparacao <- data.frame(
      Data = rep(datas_previsao, 2),
      Tipo = rep(c("Verdadeiro", "Predito"), each = h),
      Valor = c(valores_reais[[paste0("h", h)]], previsoes_multi[[paste0("h", h)]])
    )
    
    
 
    grafico_comparacao <- ggplot(dados_comparacao, aes(x = Data, y = Valor, color = Tipo, group = Tipo)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      scale_color_manual(values = c("Verdadeiro" = "#3498db", "Predito" = "#e74c3c")) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
      labs(title = paste(""),
           subtitle = paste(""),
           x = "Período", y = "Quantidade") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    print(grafico_comparacao)
    ggsave(paste0("comparacao_real_previsto_h", h, ".png"), plot = grafico_comparacao, 
           width = 10, height = 6, dpi = 300)
    
    
    dados_residuos <- data.frame(
      Data = datas_previsao,
      Residuo = residuos[[paste0("h", h)]]
    )
    
    grafico_residuos <- ggplot(dados_residuos, aes(x = Data, y = Residuo)) +
      geom_col(fill = "#f39c12", alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      labs(title = paste(""),
           x = "Data", y = "Resíduo (Verdadeiro - Predito)") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(grafico_residuos)
    ggsave(paste0("analise_residuos_h", h, "_datas.png"), plot = grafico_residuos, 
           width = 10, height = 6, dpi = 300)
  }
}


resultados_avaliacao <- data.frame()

for (h in horizontes) {
  if (h <= length(serie_temporal)) {
    mae <- mean(abs(residuos[[paste0("h", h)]]), na.rm = TRUE)
    rmse <- sqrt(mean(residuos[[paste0("h", h)]]^2, na.rm = TRUE))
    mape <- mean(abs(residuos[[paste0("h", h)]]/valores_reais[[paste0("h", h)]]), na.rm = TRUE) * 100
    
    resultados_avaliacao <- rbind(resultados_avaliacao, data.frame(
      Horizonte = h,
      MAE = mae,
      RMSE = rmse,
      MAPE = mape
    ))
  }
}

cat("\n=== AVALIAÇÃO DAS PREVISÕES COM VALIDAÇÃO WALK-FORWARD ===\n")
print(resultados_avaliacao)

print(xtable(resultados_avaliacao %>% 
         select(Horizonte, MAE, RMSE), 
       caption = "Qualidade do Ajuste de acordo com o horizonte considerado", 
       label = "tab:qualidade_ajuste_realizado"), include.rownames = F)


horizontes <- c(1, 3, 6, 12)
ultimo_valor_treino <- tail(serie_temporal, 1)


previsao_poisson_multi <- function(x_inicial, alpha, lambda, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + lambda,0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}

previsao_binomial_negativa_multi <- function(x_inicial, alpha, mu, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + mu,0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}

previsao_geometrica_multi <- function(x_inicial, alpha, mu, h) {
  previsoes <- numeric(h)
  x_atual <- x_inicial
  
  for (i in 1:h) {
    previsoes[i] <- round(alpha * x_atual + mu,0)
    x_atual <- previsoes[i]
  }
  
  return(previsoes)
}

previsoes_multi <- list()

if (melhor_modelo == "Poisson") {
  alpha <- resultados_mvc_poisson[1]
  lambda <- resultados_mvc_poisson[2]
  
  for (h in horizontes) {
    previsoes_multi[[paste0("h", h)]] <- previsao_poisson_multi(ultimo_valor_treino, alpha, lambda, h)
  }
  
} else if (melhor_modelo == "Binomial_Negativa") {
  alpha <- resultados_mvc_binneg[1]
  mu <- resultados_mvc_binneg[2]
  
  for (h in horizontes) {
    previsoes_multi[[paste0("h", h)]] <- previsao_binomial_negativa_multi(ultimo_valor_treino, alpha, mu, h)
  }
  
} else {
  alpha <- resultados_mvc_geo[1]
  mu <- resultados_mvc_geo[2]
  
  for (h in horizontes) {
    previsoes_multi[[paste0("h", h)]] <- previsao_geometrica_multi(ultimo_valor_treino, alpha, mu, h)
  }
}

modelo_poisson_estimado <- glm(serie_temporal ~ 1, family = poisson)

disp_test <- dispersiontest(modelo_poisson_estimado)

print(disp_test)

calcular_metricas <- function(previsoes, valores_reais) {
  erros <- valores_reais - previsoes
  mae <- mean(abs(erros), na.rm = TRUE)
  rmse <- sqrt(mean(erros^2, na.rm = TRUE))
  mape <- mean(abs(erros/valores_reais), na.rm = TRUE) * 100
  
  return(list(mae = mae, rmse = rmse, mape = mape))
}


resultados_avaliacao <- data.frame()

for (h in horizontes) {
  if (h <= length(serie_temporal)) {
    previsoes_h <- previsoes_multi[[paste0("h", h)]]
    valores_reais_h <- serie_temporal[1:h]
    
    metricas <- calcular_metricas(previsoes_h, valores_reais_h)
    
    resultados_avaliacao <- rbind(resultados_avaliacao, data.frame(
      Horizonte = h,
      MAE = metricas$mae,
      RMSE = metricas$rmse,
      MAPE = metricas$mape
    ))
  }
}

cat("\n=== AVALIAÇÃO DAS PREVISÕES ===\n")
print(resultados_avaliacao)
previsao_completa



serie_temporal
################################################################################
#                           PREVISÃO DE 12 PASSOS À FRENTE                     #
################################################################################

cat("\n=== PREVISÃO DE 12 PASSOS À FRENTE ===\n")
horizonte_final <- 12
ultimo_valor_treino <- tail(serie_temporal, 1)
if (melhor_modelo == "Poisson") {
  previsao_12_passos <- previsao_poisson_multi(ultimo_valor_treino, 
                                               resultados_mvc_poisson[1], 
                                               resultados_mvc_poisson[2], 
                                               horizonte_final)
} else if (melhor_modelo == "Binomial_Negativa") {
  previsao_12_passos <- previsao_binomial_negativa_multi(ultimo_valor_treino, 
                                                         resultados_mvc_binneg[1], 
                                                         resultados_mvc_binneg[2], 
                                                         horizonte_final)
} else if (melhor_modelo == "Geométrica") {
  previsao_12_passos <- previsao_geometrica_multi(ultimo_valor_treino, 
                                                  resultados_mvc_geo[1], 
                                                  resultados_mvc_geo[2], 
                                                  horizonte_final)
}

cat("Previsão dos próximos 12 passos:\n")
print(previsao_12_passos)
################################################################################
#                      GRÁFICO COM SÉRIE REAL + PREVISÃO                       #
################################################################################

ultima_data_real <- max(dados$data)
cat("Última data real:", format(ultima_data_real, "%d/%m/%Y"), "\n")
datas_previsao <- seq.Date(from = ultima_data_real %m+% months(1), 
                           by = "month", 
                           length.out = horizonte_final)

dados_completos <- data.frame(
  Data = c(dados$data, datas_previsao),
  Valor = c(serie_temporal, previsao_12_passos),
  Tipo = c(rep("Real", length(serie_temporal)), rep("Previsão", horizonte_final)),
  stringsAsFactors = FALSE
)

grafico_previsao_12 <- ggplot(dados_completos, aes(x = Data, y = Valor, color = Tipo, group = 1)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = as.numeric(ultima_data_real), 
             linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = c("Real" = "#3498db", "Previsão" = "#e74c3c")) +
  labs(title = paste(""),
       subtitle = "",
       x = "Período", 
       y = "Quantidade",
       color = "Tipo") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "red"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom")

print(grafico_previsao_12)
ggsave("previsao_12_passos_completa.png", plot = grafico_previsao_12, 
       width = 14, height = 8, dpi = 300)
################################################################################
#                         TABELA LATEX COM AS PREVISÕES                        #
################################################################################

tabela_previsoes <- data.frame(
  Data = format(datas_previsao, "%b/%Y"),
  Previsão = previsao_12_passos,
  Modelo = melhor_modelo
)

cat("\nTabela de Previsões:\n")
print(tabela_previsoes)

tabela_latex_previsoes <- xtable(tabela_previsoes %>% 
                                   select(-Modelo),
                                 caption = paste("Previsão de", horizonte_final, 
                                                 "passos à frente - Modelo", melhor_modelo),
                                 label = "tab:previsao_12_passos",
                                 digits = 0)

print(tabela_latex_previsoes, 
      file = "previsao_12_passos.tex",
      include.rownames = FALSE,
      caption.placement = "top",
      table.placement = "htbp",
      size = "small")
################################################################################
#                    GRÁFICO ZOOM NA ÁREA DE PREVISÃO                          #
################################################################################

dados_transicao <- dados_completos[(nrow(dados_completos) - horizonte_final - ultimos_reais + 1):nrow(dados_completos), ]

grafico_transicao <- ggplot(dados_transicao, aes(x = Data, y = Valor, color = Tipo, group = 1)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = as.numeric(ultima_data_real %m+% months(1)), 
             linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c("Real" = "#3498db", "Previsão" = "#e74c3c")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b%Y") +
  labs(title = paste(""),
       subtitle = paste(""),
       x = "Período", 
       y = "Quantidade",
       color = "Tipo") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(grafico_transicao)
ggsave("detalhe_transicao_previsao.png", plot = grafico_transicao, 
       width = 12, height = 6, dpi = 300)

# 10.4 ESTATÍSTICAS DESCRITIVAS DAS PREVISÕES
# ===========================================

cat("\n=== ESTATÍSTICAS DAS PREVISÕES ===\n")
cat("Média das previsões:", round(mean(previsao_12_passos), 2), "\n")
cat("Desvio padrão das previsões:", round(sd(previsao_12_passos), 2), "\n")
cat("Mínimo:", min(previsao_12_passos), "\n")
cat("Máximo:", max(previsao_12_passos), "\n")
cat("Intervalo de confiança 95%: [", 
    round(mean(previsao_12_passos) - 1.96 * sd(previsao_12_passos), 2), ", ",
    round(mean(previsao_12_passos) + 1.96 * sd(previsao_12_passos), 2), "]\n", sep = "")


################################################################################
#             ANÁLISE DO EQM DA ÚLTIMA SIMULAÇÃO DE MONTE CARLO                #
################################################################################
eqm <- data.frame(
  n = c(50, 100, 300, 500),
  
  Poisson_MQC_alpha = c(0.0019, 0.0011, 0.0005, 0.0003),
  Poisson_MQC_lambda = c(1.9345, 1.3793, 0.7112, 0.4846),
  
  Poisson_MVC_alpha = c(0.0005, 0.0002, 0.0001, 0.0000),
  Poisson_MVC_lambda = c(0.5100, 0.2789, 0.1127, 0.0688),
  
  BN_MQC_alpha = c(0.0018, 0.0011, 0.0005, 0.0003),
  BN_MQC_lambda = c(1.8906, 1.4687, 0.7348, 0.5139),
  
  BN_MVC_alpha = c(0.0010, 0.0006, 0.0002, 0.0002),
  BN_MVC_lambda = c(1.8287, 0.7740, 0.2791, 0.2555),
  
  Geo_MQC_alpha = c(0.0066, 0.0027, 0.0008, 0.0004),
  Geo_MQC_lambda = c(7.7092, 3.8210, 1.2659, 0.7191),
  
  Geo_MVC_alpha = c(0.0002, 0.0001, 0.0000, 0.0000),
  Geo_MVC_lambda = c(0.5542, 0.2771, 0.0968, 0.0573)
)

eqm_long_alpha <- eqm %>%
  select(n, contains("alpha")) %>%
  pivot_longer(cols = -n, names_to = "Metodo", values_to = "EQM") %>%
  mutate(
    Distribuicao = case_when(
      grepl("Poisson", Metodo) ~ "Poisson",
      grepl("BN", Metodo) ~ "Binomial Negativa",
      grepl("Geo", Metodo) ~ "Geométrica"
    ),
    Estimador = ifelse(grepl("MQC", Metodo), "MQC", "MVC")
  )

eqm_long_lambda <- eqm %>%
  select(n, contains("lambda")) %>%
  pivot_longer(cols = -n, names_to = "Metodo", values_to = "EQM") %>%
  mutate(
    Distribuicao = case_when(
      grepl("Poisson", Metodo) ~ "Poisson",
      grepl("BN", Metodo) ~ "Binomial Negativa",
      grepl("Geo", Metodo) ~ "Geométrica"
    ),
    Estimador = ifelse(grepl("MQC", Metodo), "MQC", "MVC")
  )


grafico_eqm_simulacao <- ggplot(eqm_long_lambda, aes(x = factor(n), y = EQM, fill = Estimador)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Distribuicao, scales = "free_y") +
  labs(title = "",
       subtitle = "",
       x = "Tamanho da amostra (n)",
       y = "EQM",
       fill = "Estimador") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")

ggsave("eqm_simulacao.png", plot = grafico_eqm_simulacao, 
       width = 12, height = 6, dpi = 300)
