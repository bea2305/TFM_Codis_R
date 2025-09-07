# Pedigree and Genetic Analysis
library(pedigreeTools)  # Pedigree analysis tools
library(nadiv)          # Tools for analyzing pedigrees and relatedness
library(dplyr)


base_path <- "C:/Users/Cristian/Desktop"
pedigreedataL <- read.csv(file.path(base_path,"PedigreeData_ind-level.csv"))
pedigreedata_ordL <- editPed(dam= pedigreedataL$dam, sire=pedigreedataL$sire, label=pedigreedataL$ID)
pedigreedata_ordL <- pedigreedata_ordL[,c(1,3,2)]
colnames(pedigreedata_ordL) <- c("ID","dam","sire")      
head(pedigreedata_ordL)

AmatL <- as.matrix(nadiv::makeA(pedigreedata_ordL))  # Crea la matriu de distàncies


lleida <- read.csv(file.path(base_path,"BreedingData_brood-level.csv"))
lleida_lh <- read.csv(file.path(base_path,"LifehistoryData_ind-level.csv"))


# Filtra los individuos que han eclosionado y guarda su año y nido de nacimiento.
lifehistdata_hatch <- lleida_lh %>% filter(stage=="hatched") %>%
  mutate(birthyear=year, # year of birth
         birthbox=nestbox.ID) %>% # nest box of birth
  dplyr::select(ID,birthyear,birthbox)


# Une los datos de puesta con los de historia vital de la hembra (fem.ID), calcula su edad en el año de puesta, y renombra fem.ID como animal
phenodataL <- left_join(lleida,lifehistdata_hatch, by=c("fem.ID"="ID"))
phenodataL <- phenodataL %>% mutate(animal=fem.ID,
                                    age=year-birthyear)

# Cuenta cuántas hembras tienen edad conocida
phenodataL %>% dplyr::filter(!is.na(age)) %>% summarise(n_distinct(animal)) # Resultado: 19 


# Crear expdata con datos de cría
expdata <- left_join(phenodataL, lleida_lh %>% filter(stage == "breeding"),
                     by = c("animal" = "ID", "year", "nestbox.ID"))


# Crea una nueva variable breed_exp que indica si la hembra es inexperta (1-2 años) o experta (>2 años). Si no hay datos, pone NA.expdata <- expdata %>%
expdata <- expdata %>%
arrange(animal, year) %>%
  mutate(
    breed_exp = case_when(
      stage == "breeding" & age %in% 1:2 ~ "unexperienced",
      stage == "breeding" & age > 2 ~ "experienced",
      TRUE ~ NA_character_))

# Corrige los años sin edad o con incertidumbre: detecta el primer año de cria para cada hembra. Cambia NA por uncertain o experienced (si no es el primero)
lifehistdata_breed <- expdata %>%
  group_by(animal) %>% 
  summarise(# Identify the first year of breeding
    first_breeding_year = case_when(stage == "breeding" & year==min(year) ~ year)) %>%
  filter(!is.na(first_breeding_year))

expdata <- left_join(expdata, lifehistdata_breed, by="animal") %>%
  mutate(
    # Assign 'uncertain' for all breeding observations in the first year of breeding
    breed_exp = case_when(is.na(breed_exp) & stage == "breeding" & year == first_breeding_year ~ "uncertain",
                          # Assign 'experienced' for all subsequent breeding years                   
                          is.na(breed_exp) & stage == "breeding" & year > first_breeding_year ~ "experienced",
                          TRUE ~ breed_exp)) %>%
  ungroup() %>%
  dplyr::select(animal,year,nestbox.ID,breed_exp) %>% distinct() %>% filter(!is.na(breed_exp))


# Añade la experiencia reproductora al dataset de puesta
phenodataL <- left_join(phenodataL,expdata,by=c("animal","year","nestbox.ID"))

# Añade información materna del pedigrí
phenodataL <- left_join(phenodataL,pedigreedata_ordL %>% dplyr::select(-sire), by=c("animal"="ID"))

# Filtra para quedarte solo con hembras del pedigrí y selecciona columnas relevantes.
phenodataL <- phenodataL %>% filter(animal %in% unique(pedigreedata_ordL$ID)) %>%
  dplyr::select(year,pop,LD,jLD,colony.ID,nestbox.ID,cl.size,fem.ID,animal,birthyear,birthbox,breed_exp,dam) %>%
  distinct()

# Resumen
nrow(phenodataL) # TOTAL: 203 obs
length(unique(phenodataL$animal, phenodataL$breed_exp)) #149 unique females
table(phenodataL$breed_exp)











# Animal Model Analysis

library(MCMCglmm)

# Es un prior no informativo. El modelo estimará la varianza residual y la genética a partir de los datos, sin imponerles casi nada.
prior1 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002)))   # No estás diciendo que la varianza sea realmente 1, solo que partes de ahí, pero con tan poca confianza (nu = 0.002) que los datos lo decidirán.

animal_model <- MCMCglmm(
  cl.size ~ breed_exp,               # Efectos fijos: experiencia reproductora
  random = ~ animal,                 # Efecto aleatorio: es el pedigrí
  pedigree = pedigreedata_ordL,      # Pedigrí ordenado
  data = phenodataL,                 # Datos fenotípicos
  family = "gaussian",               # Tipo de respuesta: continua (tamaño de puesta)
  prior = prior1,                    # Prior definido arriba
  nitt = 650000,                     # número total de iteraciones
  thin = 60,                         # cada cuántas iteraciones se guarda una muestra
  burnin = 50000)                    # primeras iteraciones que se descartan


# Calcular heredabilidad
posterior.heritability <- animal_model$VCV[,"animal"]/
  (animal_model$VCV[,"animal"] + animal_model$VCV[,"units"])   # animal_model$VCV[,"animal"]: varianza genética estimada en cada iteración. animal_model$VCV[,"units"]: varianza residual.

posterior.mode(posterior.heritability) # Valor más probable de heredabilidad de las iteraciones. Resultat: h2 = 0,6734 (valor molt alt)
HPDinterval(posterior.heritability, 0.95) # Intérvalo creíble del 95% (0,5009 a 0,7719)

# Explicación resultado: Podría sobreestimar h² si parte de la variación genética en realidad es debida a experiencias repetidas de la hembra o a años buenos/malo, ya que aquí no se está teniendo en cuenta.











#   ANIMAL MODELS WITH MAIN RELEVANT RANDOM EFFECTS

# Prior no informativo con 3 efectos aleatorios: G1 (animal, genético), G2 (year) y G3 (fem.ID, efectos maternales, efectos no genéticos individuales de las hembras)
prior4 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002),
                        G2 = list(V = 1, nu = 0.002),
                        G3 = list(V = 1, nu = 0.002)))


animalmodel_confounds <-MCMCglmm(cl.size ~ breed_exp,                           # Efectos fijos
                                 random = ~animal + year + fem.ID,              # Efectos aleatorios
                                 pedigree = pedigreedata_ordL,
                                 data = phenodataL,
                                 prior = prior4,
                                 nitt = 650000, thin = 60, burnin = 50000, verbose = TRUE) # Más iteraciones porque el modelo es más complejo

# Calcular heredabilidad
posterior.heritability <- animalmodel_confounds$VCV[,"animal"] /
  (animalmodel_confounds$VCV[,"animal"] +
     animalmodel_confounds$VCV[,"fem.ID"] +
     animalmodel_confounds$VCV[,"year"] +
     animalmodel_confounds$VCV[,"units"])

posterior.mode(posterior.heritability)           # Heredabilidad = 0,0042
HPDinterval(posterior.heritability, 0.95)        # IC = 0,0003 - 0,6615

# Explicación resultados: Hay mucha incertidumbre: el modelo no puede descartar que la heredabilidad sea moderada o incluso alta… pero también podría ser prácticamente cero.
# Efecto year: quizás los años buenos producen más huevos, y como algunas hembras solo crían en ciertos años, eso se confundía con genética.
# Efectos maternales: hay hembras que ponen más huevos, no por sus genes, sino por experiencia, condición física, edad, etc.

# Varianza explicada por el año
posterior.year <- animalmodel_confounds$VCV[,"year"] /
  (animalmodel_confounds$VCV[,"animal"] +
     animalmodel_confounds$VCV[,"fem.ID"] +
     animalmodel_confounds$VCV[,"year"] +
     animalmodel_confounds$VCV[,"units"])

posterior.mode(posterior.year)     # Resultado = 0,0026
HPDinterval(posterior.year, 0.95)  # IC = 0,0002 - 0,1790


# Varianza explicada por la hembra (no genética)
posterior.permanent_effects <- animalmodel_confounds$VCV[,"fem.ID"] /
  (animalmodel_confounds$VCV[,"animal"] +
     animalmodel_confounds$VCV[,"fem.ID"] +
     animalmodel_confounds$VCV[,"year"] +
     animalmodel_confounds$VCV[,"units"])

posterior.mode(posterior.permanent_effects)          # Resultado = 0,0046
HPDinterval(posterior.permanent_effects, 0.95)       # IC = 0,0002 - 0,6536

# Tanto year como fem.ID no explican mucho por sí solos, al menos según este modelo.
# Esto ayuda a entender por qué la heredabilidad también salió baja (~0.003): Ningún componente (genético, individual ni anual) parece explicar bien la variación en cl.size.
# Esto puede deberse a: Mucho ruido ambiental no modelado. Falta de poder estadístico (pocas hembras con datos suficientes). O que el tamaño de puesta realmente varía poco entre hembras y se ajusta mucho al ambiente inmediato.












#   ANIMAL MODELS WITH ALL RELEVANT RANDOM EFFECTS

# Igual que el anterior pero incorpora G4 (caja nido donde nace) 
prior4 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002),
                        G2 = list(V = 1, nu = 0.002),
                        G3 = list(V = 1, nu = 0.002),
                        G4 = list(V = 1, nu = 0.002)))

phenodataL_clean <- phenodataL %>% filter(!is.na(birthbox))   # Eliminar los NA de birthbox

animalmodel_earlylife_N <- MCMCglmm(
  cl.size ~ breed_exp,
  random = ~animal + year + fem.ID + birthbox,
  pedigree = pedigreedata_ordL,
  data = phenodataL_clean,
  prior = prior4,
  nitt = 650000, thin = 60, burnin = 50000,
  verbose = TRUE
)



# Calcular heredabilidad
posterior.heritability <- animalmodel_earlylife_N$VCV[,"animal"]/
  (animalmodel_earlylife_N$VCV[,"animal"] + animalmodel_earlylife_N$VCV[,"fem.ID"] + animalmodel_earlylife_N$VCV[,"year"] +
     animalmodel_earlylife_N$VCV[,"birthbox"] + animalmodel_earlylife_N$VCV[,"units"])

posterior.mode(posterior.heritability)             # Resultado = 0,0043
HPDinterval(posterior.heritability, 0.95)          # IC = 0,0002 - 0,3673


# Varianza explicada por el año
posterior.year <- animalmodel_earlylife_N$VCV[,"year"]/
  (animalmodel_earlylife_N$VCV[,"animal"] + animalmodel_earlylife_N$VCV[,"fem.ID"] + animalmodel_earlylife_N$VCV[,"year"] +
     animalmodel_earlylife_N$VCV[,"birthbox"] + animalmodel_earlylife_N$VCV[,"units"])

posterior.mode(posterior.year)                   # Resultado = 0,0030
HPDinterval(posterior.year, 0.95)                # IC = 0,0002 - 0,3215


# Varianza explicada por la hembra (no genética)
posterior.permanent_effects <- animalmodel_earlylife_N$VCV[,"fem.ID"]/
  (animalmodel_earlylife_N$VCV[,"animal"] + animalmodel_earlylife_N$VCV[,"fem.ID"] + animalmodel_earlylife_N$VCV[,"year"] +
     animalmodel_earlylife_N$VCV[,"birthbox"] + animalmodel_earlylife_N$VCV[,"units"])

posterior.mode(posterior.permanent_effects)                    # Resultado = 0,0032
HPDinterval(posterior.permanent_effects, 0.95)                 # IC = 0,0003 - 0,3869


# Varianza explicada por la caja nido en la que nació la hembra
posterior.nestbox <- animalmodel_earlylife_N$VCV[,"birthbox"]/
  (animalmodel_earlylife_N$VCV[,"animal"] + animalmodel_earlylife_N$VCV[,"fem.ID"] + animalmodel_earlylife_N$VCV[,"year"] +
     animalmodel_earlylife_N$VCV[,"birthbox"] + animalmodel_earlylife_N$VCV[,"units"])

posterior.mode(posterior.nestbox)                  # Resultado = 0,0030
HPDinterval(posterior.nestbox, 0.95)               # IC = 0,0002 - 0,4177

# Ningún componente parece explicar una gran parte de la varianza.
# Esto sugiere que el tamaño de puesta podría estar muy influenciado por factores ambientales inmediatos no incluidos en el modelo, o que hay poca variabilidad real entre hembras.










install.packages("gt")
install.packages("webshot2")  # Alternativa moderna a webshot
webshot2::install_phantomjs()  # Necesario para hacer capturas

library(gt)
library(webshot2)

# Crear la tabla como antes
resumen_modelos <- data.frame(
  Model = c("Model 1 (simple)", "Model 2 (confusió)", "Model 3 (complet)"),
  
  Heretatabilitat = c(0.6734, 0.0042, 0.0043),
  Heretabilitat_IC = c("[0.5009, 0.7719]", "[0.0003, 0.6615]", "[0.0002, 0.3673]"),
  
  Any = c(NA, 0.0026, 0.0030),
  Any_IC = c(NA, "[0.0002, 0.1790]", "[0.0002, 0.3215]"),
  
  Femella = c(NA, 0.0046, 0.0032),
  Femella_IC = c(NA, "[0.0002, 0.6536]", "[0.0003, 0.3869]"),
  
  CaixaNiu = c(NA, NA, 0.0030),
  CaixaNiu_IC = c(NA, NA, "[0.0002, 0.4177]")
)


library(gt)
tabla_gt <- resumen_modelos %>%
  gt() %>%
  tab_header(title = "Animal models") %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 4,
    drop_trailing_zeros = FALSE,
    use_seps = FALSE,
    locale = "es"
  ) %>%
  fmt_missing(columns = everything(), missing_text = "NA") %>%
  cols_align(
    align = "center",  # centra el contenido
    columns = everything()
  )

gtsave(tabla_gt, filename = "tabla_modelos.png")
































library(readr)
library(dplyr)

base_path <- "C:/Users/Cristian/Desktop"
ped <- read_csv(file.path(base_path,"PedigreeData_ind-level.csv"),
                show_col_types = FALSE)

# Neteja de NA i assegura tipus
ped <- ped %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(c(dam, sire), ~na_if(., ""))) %>%
  mutate(across(c(dam, sire), ~ifelse(. %in% c("NA","NaN","nan","None"), NA, .)))

id_set <- ped$ID
dam_map <- setNames(ped$dam, ped$ID)
sire_map <- setNames(ped$sire, ped$ID)

gen_cache <- new.env(parent = emptyenv())

get_gen <- function(id, stack = character()) {
  if (!is.null(gen_cache[[id]])) return(gen_cache[[id]])
  if (id %in% stack) { # prevenció de cicles
    gen_cache[[id]] <- 0
    return(0)
  }
  d <- dam_map[[id]]
  s <- sire_map[[id]]
  # Pare/mare fundadors si no consten a la taula
  pg <- integer(0)
  for (p in c(d, s)) {
    if (is.na(p) || !(p %in% id_set)) {
      pg <- c(pg, 0)
    } else {
      pg <- c(pg, get_gen(p, c(stack, id)))
    }
  }
  g <- if (length(pg) == 0) 0 else 1 + max(pg)
  gen_cache[[id]] <- g
  g
}

gens <- vapply(id_set, get_gen, FUN.VALUE = numeric(1))
ped_gen <- ped %>% mutate(generation = gens)

# Resultats
max_gen <- max(ped_gen$generation)
dist <- sort(table(ped_gen$generation))

max_gen
dist
head(ped_gen)

