library(purrr)
library(broom)
library(broom.mixed)
library(lmerTest) 
library(dplyr)
library(gridExtra)
library(grid)



prod <- read.csv("/Users/Cristian/Desktop/Productivity_FV.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
prod <- subset(prod, prod$year>2014)


# Selección por nombres (incluye n_fledge)
fitness <- prod %>%
  select(year, siteID, nestboxID, jLD, n_visits, clutch_size, n_fledge) %>%
  mutate(
    jLD = as.numeric(jLD),
    n_visits = as.numeric(n_visits),
    clutch_size = as.numeric(clutch_size),
    n_fledge = as.numeric(n_fledge)
  ) %>%
  tidyr::drop_na(clutch_size, n_fledge) %>%   # quita NA en variables clave
  mutate(fyear = factor(year))


# Calcular fitness relativa anual, dividiendo los fledging por la media de cada año

prod_filtrado <- fitness %>%
  group_by(year) %>%
  mutate(
    fitness_rel = n_fledge / mean(n_fledge, na.rm = TRUE)  # compute relative fitness
  ) %>%
  ungroup() %>%
  mutate(
    clutch_size_std = scale(clutch_size),  # standardize clutch size
    jLD_std = scale(jLD)                   # standardize jLD
  ) # Volver a dejar los datos como estaban, sin agrupar, para los siguientes cálculos

prod_filtrado <- data.frame(prod_filtrado)




resultats_models <- prod_filtrado %>%
  group_split(year) %>%
  map_df(function(df) {
    yr <- unique(df$year)
    
    # Model 1: S_dif
    m1 <- lm(fitness_rel ~ clutch_size_std, data = df)
    coef1 <- tidy(m1) %>%
      filter(term == "clutch_size_std") %>%
      mutate(tipus = "S_dif", model = "model1", year = yr)
    
    # Model 2: S_dif + C_dif
    m2 <- lm(fitness_rel ~ clutch_size_std + I(clutch_size_std^2), data = df)
    coef2 <- tidy(m2) %>%
      filter(term %in% c("clutch_size_std", "I(clutch_size_std^2)")) %>%
      mutate(tipus = ifelse(term == "clutch_size_std", "S_dif", "C_dif"),
             model = "model2", year = yr)
    
    # Model 3: S_grad + C_grad
    m3 <- lmer(fitness_rel ~ clutch_size_std + I(clutch_size_std^2) + jLD_std + (1 | siteID),
               data = df, REML = FALSE)
    coef3 <- tidy(m3, effects = "fixed") %>%
      filter(term %in% c("clutch_size_std", "I(clutch_size_std^2)")) %>%
      mutate(tipus = ifelse(term == "clutch_size_std", "S_grad", "C_grad"),
             model = "model3", year = yr)
    
    bind_rows(coef1, coef2, coef3)
  })

taula_model1 <- resultats_models %>%
  filter(model == "model1") %>%
  select(year, term, estimate, std.error, p.value, tipus) %>%
  mutate(term = recode(term,
                       "clutch_size_std" = "mida_posta"),
         across(where(is.numeric), \(x) round(x, 3)))

library(clipr)

write_clip(taula_model1)


taula_model2 <- resultats_models %>%
  filter(model == "model2") %>%
  select(year, term, estimate, std.error, p.value, tipus) %>%
  mutate(term = recode(term,
                       "clutch_size_std" = "mida_posta",
                       "I(clutch_size_std^2)" = "mida_posta^2"),
         tipus = recode(tipus,
                        "C_dif" = "Q_dif"),
         across(where(is.numeric), \(x) round(x, 3)))

write_clip(taula_model2)

taula_model3 <- resultats_models %>%
  filter(model == "model3") %>%
  select(year, term, estimate, std.error, p.value, tipus) %>%
  mutate(term = recode(term,
                       "clutch_size_std" = "mida_posta",
                       "I(clutch_size_std^2)" = "mida_posta^2"),
         tipus = recode(tipus,
                        "C_grad" = "Q_grad"),
         across(where(is.numeric), \(x) round(x, 3)))

write_clip(taula_model3)


# Guardar taula del Model 1 – Diferencial direccional (S_dif)
png("taula_model1_dif_direccional.png", width = 1000, height = 600)
grid.draw(tableGrob(taula_model1))
dev.off()

# Guardar taula del Model 2 – Diferencials quadràtics (S_dif i C_dif)
png("taula_model2_dif_quadratic.png", width = 1000, height = 600)
grid.draw(tableGrob(taula_model2))
dev.off()

# Guardar taula del Model 3 – Gradients de selecció (S_grad i C_grad)
png("taula_model3_gradients.png", width = 1000, height = 600)
grid.draw(tableGrob(taula_model3))
dev.off()









library(ggplot2)

graf1 <- ggplot(taula_model1, aes(x = year, y = estimate, color = tipus, shape = tipus)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = 2015:2025) +
  labs(
    x = "Any",
    y = "Diferencial de selecció ± SE",
  ) +
  theme_minimal(base_size = 14)
dev.off()  
print(graf1)



ggplot(taula_model2, aes(x = year, y = estimate, color = tipus, shape = tipus)) +
  geom_point(size = 3, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = 2015:2025) +
  labs(
    x = "Any",
    y = "Diferencial de selecció ± SE",
  ) +
  theme_minimal(base_size = 14)


ggplot(taula_model3, aes(x = year, y = estimate, color = tipus, shape = tipus)) +
  geom_point(size = 3, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = 2015:2025) +
  labs(
    x = "Any",
    y = "Gradient de selecció ± SE",
  ) +
  theme_minimal(base_size = 14)




ggplot(prod_filtrado, aes(x = clutch_size, y = n_fledge)) +
  geom_smooth(method = "loess", se = TRUE, color = "darkgreen", linewidth = 1) +
  facet_wrap(~ year, scales = "free") +
  scale_x_continuous(
    breaks = function(x) pretty(x)[pretty(x) %% 1 == 0]  # Solo enteros
  ) +
  labs(
    x = "Mida de posta",
    y = "Nombre de polls volanders"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12)
  ) 
