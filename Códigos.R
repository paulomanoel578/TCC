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
#                          SIMULAÇÕES UM PASSO A FRENTE                        #
################################################################################
resultado_previsao <- simulacao_monte_carlo_previsao_poisson(
  n_simulacoes = 5000, 
  n_treino = 40, 
  n_teste = 10, 
  alpha_verdadeiro = 0.1, 
  lambda_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_poisson(
  n_simulacoes = 5000, 
  n_treino = 80, 
  n_teste = 20, 
  alpha_verdadeiro = 0.1, 
  lambda_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_poisson(
  n_simulacoes = 5000, 
  n_treino = 240, 
  n_teste = 60, 
  alpha_verdadeiro = 0.1, 
  lambda_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_poisson(
  n_simulacoes = 5000, 
  n_treino = 400, 
  n_teste = 100, 
  alpha_verdadeiro = 0.1, 
  lambda_verdadeiro = 5
)


print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

################################################################################
#                             BINOMINAL NEGATIVA                               #
################################################################################

resultado_previsao <- simulacao_monte_carlo_previsao_bin_negativa(
  n_simulacoes = 5000, 
  n_treino = 40, 
  n_teste = 10, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5, 
  sigma_verdadeiro = 6
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_bin_negativa(
  n_simulacoes = 5000, 
  n_treino = 80, 
  n_teste = 20, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5, 
  sigma_verdadeiro = 6
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_bin_negativa(
  n_simulacoes = 5000, 
  n_treino = 240, 
  n_teste = 60, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5, 
  sigma_verdadeiro = 6
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_bin_negativa(
  n_simulacoes = 5000, 
  n_treino = 400, 
  n_teste = 100, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5, 
  sigma_verdadeiro = 6
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

################################################################################
#                             GEOMETRICA                                       #
################################################################################

resultado_previsao <- simulacao_monte_carlo_previsao_geometrica(
  n_simulacoes = 5000, 
  n_treino = 40, 
  n_teste = 10, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_geometrica(
  n_simulacoes = 5000, 
  n_treino = 80, 
  n_teste = 20, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_geometrica(
  n_simulacoes = 5000, 
  n_treino = 240, 
  n_teste = 60, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

resultado_previsao <- simulacao_monte_carlo_previsao_geometrica(
  n_simulacoes = 5000, 
  n_treino = 400, 
  n_teste = 100, 
  alpha_verdadeiro = 0.1, 
  mu_verdadeiro = 5
)

print(paste("MAE médio CLS:", resultado_previsao$mae_media_CLS))
print(paste("RMSE médio CLS:", resultado_previsao$rmse_media_CLS))
print(paste("MAE médio CML:", resultado_previsao$mae_media_CML))
print(paste("RMSE médio CML:", resultado_previsao$rmse_media_CML))

