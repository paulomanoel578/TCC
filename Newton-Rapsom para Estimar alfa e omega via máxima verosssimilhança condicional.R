rm(list=ls(all=T))
gc()


# Simulação de dados do modelo Poisson INAR(1)
set.seed(123)

n <- 1000  # Tamanho da amostra
alpha <- 0.15  # Parâmetro verdadeiro de autocorrelação
lambda <- 20   # Parâmetro verdadeiro de Poisson

Y <- numeric(n)
Y[1] <- rpois(1, lambda / (1 - alpha))

for (t in 2:n) {
  Y[t] <- rbinom(1, Y[t - 1], alpha) + rpois(1, lambda)
}

# Visualizando os primeiros valores
Y[1:10]

# Função de verossimilhança condicional
log_likelihood <- function(params, Y) {
  alpha <- params[1]
  lambda <- params[2]
  n <- length(Y)
  
  logL <- 0
  for (t in 2:n) {
    mu <- alpha * Y[t - 1] + lambda
    logL <- logL + dpois(Y[t], mu, log = TRUE)
  }
  return(-logL)  # Retorna o negativo da log-verossimilhança para minimização
}

# Método de Newton-Raphson para otimização
newton_raphson <- function(Y, init_params, tol = 1e-8, max_iter = 100) {
  params <- init_params
  n_iter <- 0
  
  while (n_iter < max_iter) {
    n_iter <- n_iter + 1
    
    # Calculando a log-verossimilhança e o gradiente (derivadas)
    ll <- log_likelihood(params, Y)
    grad <- numDeriv::grad(log_likelihood, params, Y = Y)
    hess <- numDeriv::hessian(log_likelihood, params, Y = Y)
    
    # Atualizando os parâmetros
    delta <- solve(hess, grad)
    params <- params - delta
    
    # Verificando a convergência
    if (sqrt(sum(delta^2)) < tol) {
      cat("Convergência atingida na iteração", n_iter, "\n")
      break
    }
  }
  
  list(params = params, iterations = n_iter, tol = sqrt(sum(delta^2)))
}

# Estimativa inicial dos parâmetros
init_params <- c(0.1, 20)  # Valores iniciais para alfa e lambda

# Rodando o método de Newton-Raphson
result <- newton_raphson(Y, init_params)

# Mostrando os resultados
cat("Estimativas de Alfa e Lambda:\n")
cat("Alfa:", result$params[1], "\n")
cat("Lambda:", result$params[2], "\n")
cat("Tolerância", result$tol)
