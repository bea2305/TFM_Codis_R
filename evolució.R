library(dplyr)
library(ggplot2)
library(brms)


prod <- read.csv("/Users/Cristian/Desktop/Productivity_FV.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
prod <- subset(prod, prod$year>2000)
prod[] <- lapply(prod, function(x) if(is.factor(x)) factor(x) else x)
str(prod)


prod <- prod %>%
  filter(
    !year %in% c(1998, 1999)
  )

table(prod$year)

prod$clutch_size <- as.numeric(prod$clutch_size)



# Prepare data and plot clutch size over the years

summary_data <- prod %>%
  group_by(year) %>%
  summarise(
    mean_clutch = mean(clutch_size, na.rm = TRUE),
    se_clutch = sd(clutch_size, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

ggplot(summary_data, aes(x = factor(year), y = mean_clutch)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_clutch - se_clutch,
                    ymax = mean_clutch + se_clutch), width = 0.2) +
  geom_text(aes(label = n, y = mean_clutch + se_clutch + 0.05), size = 3) +  # add sample size
  labs(x = "Any", y = "Mida de posta mitjana ± SE") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Formal analysis

m1 <- brm(clutch_size ~ n_visits + year + (1|siteID/nestboxID) + (1|year),
          data = prod,
          family=gaussian(),
          control = list(adapt_delta = 0.9999, max_treedepth = 12),
          cores = 6,  # cores cannot be higher than chains
          warmup=1000, iter=2000, chains=2)

summary(m2)

library(dplyr)
library(gridExtra)
library(grid)

# Extraer resumen como data.frame
tab <- as.data.frame(round(summary(m1)$fixed, 3))

# Cambiar nombres de filas (predictores)
tab <- tab %>%
  tibble::rownames_to_column("Predictor") %>%   # pasa rownames a columna
  mutate(Predictor = recode(Predictor,
                            "n_visits" = "n_visites",
                            "year" = "any")) %>%
  column_to_rownames("Predictor")   # volver a ponerlo como rownames si quieres

# Cambiar nombres de columnas
colnames(tab) <- recode(colnames(tab),
                        "Estimate" = "Valor estimat",
                        "Est.Error" = "Error estàndard")

# Pasar a tabla y guardar
tabla <- tableGrob(tab)
png("tabla_coeficients.png", width = 1000, height = 600)
grid.draw(tabla)
dev.off()




