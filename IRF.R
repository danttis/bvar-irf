
# Pacotes
library(vars)
library(BVAR)
library(BVARverse)
library(dplyr)
library(ggplot2)


# Dados
data <- fred_qd[, c("CPIAUCSL", "UNRATE")]
data <- fred_transform(data, codes = c(5, 5), lag = 4)


# VAR clássico (vars)
var_model <- VAR(y = data, p = 1, type = "const")

var_irf <- vars::irf(
  var_model,
  impulse = "CPIAUCSL",
  response = "UNRATE",
  n.ahead = 12,
  boot = FALSE
)

# Transformar VAR clássico em data.frame para ggplot
var_irf_df <- data.frame(
  time = 1:nrow(var_irf$irf$CPIAUCSL),
  irf  = var_irf$irf$CPIAUCSL[, "UNRATE"]
)


# BVAR
bvar_model <- bvar(
  data,
  lags = 1,
  n_draw = 6000L,
  n_burn = 100L,
  verbose = FALSE
)

bvar_irf <- irf(
  bvar_model,
  horizon = 12,
  identification = TRUE,  
  fevd = FALSE
)

irf_bvar <- augment(bvar_irf)

bvar_plot <- irf_bvar %>%
  filter(impulse == "CPIAUCSL", response == "UNRATE") %>%
  select(time, q16, q50, q84) %>%
  mutate(model = "BVAR")


# Gráfico comparativo VAR vs BVAR
ggplot() +
  # Linha BVAR
  geom_line(
    data = bvar_plot,
    aes(x = time, y = q50, color = "BVAR"),
    size = 1
  ) +
  # Faixa de credibilidade BVAR
  geom_ribbon(
    data = bvar_plot,
    aes(x = time, ymin = q16, ymax = q84),
    fill = "blue",
    alpha = 0.2
  ) +
  # Linha VAR clássico
  geom_line(
    data = var_irf_df,
    aes(x = time, y = irf, color = "VAR"),
    size = 1,
    linetype = "dashed"
  ) +
  scale_color_manual(values = c("BVAR" = "blue", "VAR" = "red")) +
  labs(
    title = "Resposta de UNRATE a choque em CPIAUCSL",
    x = "Horizonte",
    y = "IRF",
    color = "Modelo"
  ) +
  theme_minimal()
