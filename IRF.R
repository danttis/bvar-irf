# Pacotes necessários
library(vars)
library(MASS)  # Para mvrnorm

set.seed(123)

# --- Dados simulados ---
K <- 2    # número de variáveis
T <- 100  # observações
H <- 10   # horizonte IRF
p <- 1    # lag

# Simulando um VAR(1)
A_true <- matrix(c(0.5, 0.2,
                   0.1, 0.4), nrow=K, byrow=TRUE)
Sigma_true <- matrix(c(1,0.3,
                       0.3,1), nrow=K)

Y <- matrix(0, nrow=T, ncol=K)
eps <- mvrnorm(T, mu=rep(0,K), Sigma=Sigma_true)

for(t in 2:T){
  Y[t,] <- A_true %*% Y[t-1,] + eps[t,]
}

colnames(Y) <- paste0("Var", 1:K)  # nomes para VAR

# --- 1) IRF com VAR clássico ---
var_model <- VAR(Y, p=p, type="const")
irf_var <- vars::irf(var_model, n.ahead=H, ortho=TRUE, boot=100)

# --- 2) BVAR simplificado (R puro) ---
S <- 500000
posterior_coefs <- lapply(1:S, function(s) A_true + matrix(rnorm(K^2,0,0.05), K, K))
posterior_sigma <- lapply(1:S, function(s) Sigma_true + matrix(rnorm(K^2,0,0.02), K, K))

# Funções BVAR
compute_ma_matrices <- function(A_list, p, K, H){
  comp_mat <- matrix(0, nrow=K*p, ncol=K*p)
  comp_mat[1:K, ] <- A_list
  if(p>1) comp_mat[(K+1):(K*p), 1:(K*(p-1))] <- diag(K*(p-1))
  J <- cbind(diag(K), matrix(0, K, K*(p-1)))
  Psi <- array(0, dim=c(K,K,H+1))
  Psi[,,1] <- diag(K)
  for(h in 1:H){
    comp_mat <- comp_mat %*% comp_mat
    Psi[,,h+1] <- J %*% comp_mat %*% t(J)
  }
  Psi
}

compute_irf_single <- function(A_list, Sigma, H){
  K <- nrow(Sigma)
  Psi <- compute_ma_matrices(A_list, 1, K, H)
  B <- t(chol(Sigma))
  IRF <- array(0, dim=c(K,K,H+1))
  for(h in 0:H){
    IRF[,,h+1] <- Psi[,,h+1] %*% B
  }
  IRF
}

compute_bvar_irf <- function(posterior_coefs, posterior_sigma, H){
  S <- length(posterior_coefs)
  K <- nrow(posterior_sigma[[1]])
  IRFs <- array(NA, dim=c(S,K,K,H+1))
  for(s in 1:S){
    IRFs[s,,,] <- compute_irf_single(posterior_coefs[[s]], posterior_sigma[[s]], H)
  }
  list(
    mean = apply(IRFs, c(2,3,4), mean),
    q05  = apply(IRFs, c(2,3,4), quantile, 0.05),
    q95  = apply(IRFs, c(2,3,4), quantile, 0.95)
  )
}

bvar_irfs <- compute_bvar_irf(posterior_coefs, posterior_sigma, H)

# --- 3) Plot completo BVAR vs VAR (VAR em preto) ---
par(mfrow=c(K,K), mar=c(4,4,2,1))  # grid KxK

cores <- c("blue","red")  # cores para choques do BVAR

for(i in 1:K){         # variável resposta
  for(j in 1:K){       # choque
    # BVAR
    bvar_mean <- bvar_irfs$mean[i,j,]
    bvar_q05  <- bvar_irfs$q05[i,j,]
    bvar_q95  <- bvar_irfs$q95[i,j,]
    
    # VAR
    var_mean  <- irf_var$irf[[i]][,j]
    var_lower <- irf_var$Lower[[i]][,j]
    var_upper <- irf_var$Upper[[i]][,j]
    
    # limites y
    ymin <- min(c(bvar_q05, var_lower))
    ymax <- max(c(bvar_q95, var_upper))
    
    # Plot BVAR com sombra
    plot(0:H, bvar_mean, type='l', col=cores[j], lwd=2,
         ylim=c(ymin, ymax),
         xlab="Horizon", ylab="Response",
         main=paste0("Var", i, " <- Shock ", j))
    polygon(c(0:H, rev(0:H)),
            c(bvar_q05, rev(bvar_q95)),
            col=adjustcolor(cores[j], alpha.f=0.2), border=NA)
    lines(0:H, bvar_mean, col=cores[j], lwd=2)
    
    # Sobrepor VAR em preto
    lines(0:H, var_mean, col="black", lwd=2)
    lines(0:H, var_lower, col="black", lty=3)
    lines(0:H, var_upper, col="black", lty=3)
    
    # Legenda
    legend("topright",
           legend=c("BVAR mean","BVAR CI","VAR mean","VAR CI"),
           col=c(cores[j], cores[j], "black","black"),
           lty=c(1,1,1,3), lwd=2, cex=0.7, bg="white")
  }
}

