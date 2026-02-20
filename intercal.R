#' ---
#' title: "Analyses of detector richness variance and between-detector beta diversity"
#' author: "Luciano L.M. De Benedictis and Marco Cervellini"
#' output: 
#'  pdf_document:
#'   latex_engine: xelatex
#' mainfont: DejaVu Sans
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE
)

#' # Setup
#' Set up libraries, options, import and process data

# setup --------------------------------------------------------------
library(tidyverse)
library(readxl)
library(patchwork)
library(brms)
library(bayesplot)
library(distributional)
library(tidybayes)
library(ggdist)
library(vegan)
library(patchwork)

options(brms.file_refit = "on_change")
options(mc.cores = parallel::detectCores())
theme_set(theme_minimal())

if (!dir.exists("fits")) dir.create("fits")

set.seed(932)

## Plotting functions ------------------------------------------------------

post_prior_beta <- function(x){
  x |>
    as_draws_df("b_exercise2") |> 
    ggplot(aes(x = b_exercise2)) +
    stat_slab(aes(xdist = dist_normal(0, 1)),
              fill = "#00BA38", alpha = 0.5, inherit.aes = F,
              data = tibble()) +
    geom_density(fill = "#F8766D", color = NA, alpha = 0.5)
}

post_slab <- function(data, x, line){
  data |> 
    ggplot(aes(x = {{x}})) +
    stat_slab(aes(fill_ramp = after_stat(level)),
            .width = c(.50, .80, .95, 1),
            slab_color = "navy",
            fill = "skyblue") +
    scale_fill_ramp_discrete() +
    guides(fill_ramp = "none") +
    geom_vline(xintercept = line, linetype = 2) +
    stat_pointinterval(aes(interval_color = after_stat(level)),
                       point_interval = "median_qi",
                       side = "top",
                       .width = c(.50, .80, .95),
                       interval_size_domain = c(1.5, 5)) +
    scale_color_discrete(palette = c("grey60", "grey40", "black"),
                         aesthetics = "interval_color",
                         guide = NULL) +
    theme(axis.text.y = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

plot_post <- function(data, line = 0, fit = F){
  if(fit){
    data <- data |> 
      as_draws_df() |> 
      select(value = b_exercise2) 
  }
  post_slab(data, value, line)
  }

plot_comp <- function(x, comparison){
  x |> 
    epred_draws(newdata = tibble(exercise = as.character(1:2)),
                re_formula = NA) |>
    ungroup() |> 
    select(.draw, exercise, .epred) |>
    pivot_wider(names_from = exercise, values_from = .epred, 
                names_prefix = "ex") |>
    mutate(contrast = ex2 - ex1,
         lift = ex2 / ex1 - 1) |> 
    post_slab({{comparison}}, line = 0)
}

plot_pred <- function(x){
  tibble(exercise = as.character(1:2),
         SU = NA) |> 
    add_epred_draws(x, re_formula = NA) |> 
    ggplot(aes(x = .epred, y = exercise)) +
    stat_dotsinterval(aes(interval_color = after_stat(level)),
      quantiles = 100,
      .width = c(0.5, 0.8, 0.95),
      interval_size_domain = c(1, 6),
      slab_fill = "skyblue",
      slab_color = NA) +
    scale_color_discrete(palette = c("grey60", "grey40", "black"),
                         aesthetics = "interval_color",
                         guide = NULL) +
    scale_y_discrete(limits = rev, 
                     expand = expansion(mult = c(.15, 0), add = c(0, -1)))
}

#+ ## Importing data

# data import -------------------------------------------------------------

## 2023 --------------------------------------------------------------------

species_2023 <- read_xls("data/2023/23pre_post_training_sp.xls") |> 
  select(event = 'Numero Plot', SU = 'Subplot', Layer = 'Codice Strato',
         species = 'Genere e Specie Pignatti', cover = 'Copertura %')

#recode event into detector and exercise
lookup <- tibble(event = 32:39,
                 exercise = rep(1:2, each = 4),
                 detector = rep(LETTERS[1:4], times = 2))

species_2023 <- species_2023 |> 
  left_join(lookup) |> 
  select(!event) |> 
  rename(Detector = detector)

#convert 2023 species layers

#' layers coded as:
#' 
#' 1. tree (woody + pertinent lianas and climbers ) h > 5 m
#' 2. shrub (woody + pertinent lianas and climbers ) 0.5 < h ≤ 5
#' 3. herb(however herbaceous and ferns; woody only if h ≤ 0.5 m) h ≤ 0.5 m
#' 4. moss (terricolous bryophytes and t-lichens) on mineral, organic soil

species_2023 <- species_2023 |> 
  filter(Layer !=4) |>  # 2025 doesn't have moss
  mutate(Layer = factor(Layer, labels = c("T", "S", "H")))

## 2025 --------------------------------------------------------------------

### 1° exercise ------------------------------------------------------------

species_2025_bf <- list.files("data/2025/25_02", full.names = T) |>
  map(\(x) read_excel(x, sheet = "Species")) |> 
  list_rbind() |> 
  select(!Notes) |> 
  mutate(SU = as_factor(SU)) |> 
  rename(species = "species (Genus + species)",
         cover = "Cover(perc.)") |> 
  mutate(cover = as.numeric(gsub(",", ".", cover))) |> 
  filter(Layer %in% c("T", "S", "H")) |> #layer codes
  separate_wider_delim(Plot.code, delim = ".", names = c("Detector", "Anno")) |> 
  select(!Anno) 

### 2° exercise ------------------------------------------------------------

second_2025 <- list.files("data/2025/27_02", full.names = T)
species_sheets <- second_2025[str_detect(second_2025, regex("species", ignore_case = T))]

species_2025_af <- lapply(species_sheets, function(x) read_csv(x, col_names = F, skip = 1)) |>
  bind_rows() |> 
  select(id = X1, SU = X2, Layer=X3, species=X4, cover=X5) |> 
  separate_wider_delim(id, delim = ".", names = c("Detector", "Anno")) |> 
  select(!Anno)|> 
  mutate(SU = as_factor(SU)) 

### join the exercises -----------------------------------------------------

species_2025 <- bind_rows(species_2025_bf, species_2025_af, .id = "exercise") |> 
  filter(cover>0) #f and c inserted "0" as cover obtaining all the species...

## uniform column types ---------------------------------------------------

species_2023 <- species_2023 |> 
  mutate(SU = factor(SU),
         exercise = as.character(exercise),
         Layer = as.character(Layer))

## check species names -----------------------------------------------------

# remove leading and trailing whitespaces
species_2023$species <- trimws(species_2023$species)

species_2025$species <- trimws(species_2025$species)

# 2023
unique(species_2023$species) |> sort() #sort and check all the species surveyed

# unify and correct differences among species names

species_2023[species_2023$species == "Arisarum cfr. vulgare", "species"] = "Arisarum vulgare"

species_2023[species_2023$species == "Rosa cfr. sempervirens", "species"] = "Rosa sempervirens"

# 2025
unique(species_2025$species) |> sort()

# correct species names

species_2025 <- species_2025 |> 
  # collapse cfr. to most likely species
  mutate(species = str_remove(species, "cfr\\.|cfr")) |> 
  # remove subsp. and var.
  mutate(species = str_remove(species, "(var|subsp)\\..*$")) |> 
  # change double space to single
  mutate(species = str_squish(species))

### fix typos ---------------------------------------------------------------

#' here the goal is to fix only those typos that result in the same species being recorded with different names. Typos and notes that don't lead to this are kept

species_2025$species |> unique() |> sort()

species_2025 <- species_2025 |> 
  mutate(species = str_replace_all(species,
                               c("Anemone" = "Anemonoides",
                                 "Arcticum" = "Arctium",
                                 "nemorosus" = "nemorosum",
                                 "Crategus" = "Crataegus",
                                 "Phyllirea" = "Phillyrea",
                                 "domesticus" = "domestica",
                                 "Surbus" = "Sorbus",
                                 "Symphytum sp$" = "Symphytum sp.",
                                 "Gallium" = "Galium",
                                 "lathyrus" = "Lathyrus",
                                 "Rubus sp$" = "Rubus sp.",
                                 "hulmifolius" = "ulmifolius",
                                 "Smilas" = "Smilax")))
  
species_2025$species |> unique() |> sort()

# calculate beta ---------------------------------------------------------

# function from long df to input for vegdist
vegan_cook <- function(x){
  mat <- x |> 
    mutate(cover = log1p(cover)) |> 
    pivot_wider(names_from = species, values_from = cover,
                #some detectors reported the same species twice, use the first one
                values_fn = ~ .x[[1]], values_fill = 0) |> 
    column_to_rownames("Detector") |> 
    as.matrix()
}

# function to map vegdist and wrangle output

tidyvegdist <- function(x, method, binary){
  x |> 
    map(\(x) vegdist(x, method = method, binary = binary) |> 
          as.matrix() |> 
          `diag<-`(NA) |> 
          as_tibble() |>
          pivot_longer(everything(), names_to = "detector") |> 
          filter(!(is.na(value)))) |> 
    bind_rows(.id = "ex.SU") |> 
    separate_wider_delim(ex.SU, delim = ".", names = c("exercise", "SU")) |> 
    rename(dissimilarity = value)
}

#dataframes for beta diversity
#only herb layer, group for group_split later
vegans_2025 <- species_2025 |> 
  filter(Layer == "H") |> 
  select(!c(Layer)) |> 
  group_by(exercise, SU)

#detector m didn't sample species before training in SU 7, remove also from after
vegans_2025 <- vegans_2025 |> 
  filter(!(Detector == "m" & SU == "7"))

vegans_2023 <- species_2023 |> 
  filter(Layer == "H") |> 
  select(!c(Layer)) |> 
  group_by(exercise, SU)

#list names
keys_2025 <- vegans_2025 |> 
  group_keys()

keys_2025 <- paste(keys_2025$exercise, keys_2025$SU, sep = ".")

keys_2023 <- vegans_2023 |> 
  group_keys()

keys_2023 <- paste(keys_2023$exercise, keys_2023$SU, sep = ".")

#make list of matrices for each exercise and SU
vegans_2025 <- vegans_2025 |> 
  group_split(.keep = F) |> 
  map(vegan_cook)

names(vegans_2025) <- keys_2025

vegans_2023 <- vegans_2023 |> 
  group_split(.keep = F) |> 
  map(vegan_cook)

names(vegans_2023) <- keys_2023

rm(keys_2023, keys_2025)

## Jaccard -------------------------------------------------------------

jaccard_2025 <- vegans_2025 |> 
  tidyvegdist(method = "jaccard", binary = T)

jaccard_2023 <- vegans_2023 |> 
  tidyvegdist(method = "jaccard", binary = T)

## Bray ----------------------------------------------------------------

bray_2025 <- vegans_2025 |> 
  tidyvegdist(method = "bray", binary = F)

bray_2023 <- vegans_2023 |> 
  tidyvegdist(method = "bray", binary = F)

## Euclidean -----------------------------------------------------------

euc_2025 <- vegans_2025 |> 
  tidyvegdist(method = "euclidean", binary = F) |> 
  mutate(dissimilarity = dissimilarity/max(dissimilarity))

euc_2023 <- vegans_2023 |> 
  tidyvegdist(method = "euclidean", binary = F) |> 
  mutate(dissimilarity = dissimilarity/max(dissimilarity))

## nudge 0 and 1s for beta regression ----------------------------------

jaccard_2025 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

jaccard_2023 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

bray_2025 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

bray_2023 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

euc_2025 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

euc_2023 |> 
  filter(dissimilarity == 0 | dissimilarity == 1)

euc_2025 <- euc_2025 |> 
  mutate(dissimilarity = pmin(1 - 1e-06, pmax(1e-06, dissimilarity)))

euc_2023 <- euc_2023 |> 
  mutate(dissimilarity = pmin(1 - 1e-06, pmax(1e-06, dissimilarity)))

#' # Models based on richness

# richness models ------------------------------------------------------------

richness_2023 <- species_2023 |> 
  group_by(Detector, SU, exercise) |> 
  distinct(species) |> 
  count() |> 
  mutate(Detector = as.factor(paste(Detector, exercise, sep = ".")), 
         exercise = as.factor(exercise))

richness_2025 <- species_2025 |> 
  group_by(Detector, SU, exercise) |> 
  distinct(species) |> 
  count() |> 
  mutate(Detector = as.factor(paste(Detector, exercise, sep = ".")), 
         exercise = as.factor(exercise))


# formula

rich_formula <- brmsformula(n ~ exercise + SU + (1 | gr(Detector, by = exercise)))

default_prior(rich_formula, data = richness_2025)

## priors------------------------------------------------------------------

#default prior
ggplot()+
  stat_halfeye(aes(xdist = dist_truncated(dist_student_t(3, 0, 4.4), lower = 0)),
               fill = "thistle", .width = c(0.66, 0.95))+
  scale_y_continuous(NULL, breaks = NULL)+
  ggtitle(expression("Default prior for"~sigma))

#adjusted priors
ggplot()+
  stat_halfeye(aes(xdist = dist_truncated(dist_normal(0, 5), lower = 0)),
               fill = "cornflowerblue", .width = c(0.66, 0.95))+
  scale_y_continuous(NULL, breaks = NULL)+
  ggtitle(expression("Adjusted prior for"~sigma))

ggplot()+
  stat_halfeye(aes(xdist = dist_normal(0, 5)),
               fill = "chocolate1", .width = c(0.66, 0.95))+
  scale_y_continuous(NULL, breaks = NULL)+
  ggtitle(expression("Adjusted prior for"~beta))

prior_rich <- c(prior(normal(0, 5), class = "b"), #effect of exercise and SU
                prior(normal(0, 5), class = "sd", lb = 0)) #multilevel SD


## 2025 ------------------------------------------------------------------

brm_rich_2025 <- brm(formula = rich_formula,
               data = richness_2025,
               prior = prior_rich,
               file = "fits/rich_2025",
               sample_prior = "yes",
               seed = 2958)

#check prior used
prior_summary(brm_rich_2025)

#prior predictive check for sd

prior_draws(brm_rich_2025) |> 
  select(sd_Detector) |>
  mcmc_areas(border_size = 1, prob = 0.95, point_est = "median")

#posterior summary
summary(brm_rich_2025)

brm_rich_2025 |> 
  mcmc_intervals(pars = vars(b_exercise2:sigma))

#calculate contrast
posterior_2025 <- as_draws_df(brm_rich_2025) |> 
  select(sd_before = 'sd_Detector__Intercept:exercise1',
         sd_after = 'sd_Detector__Intercept:exercise2') |> 
  mutate(contrast = sd_after - sd_before,
         ratio = sd_after^2/sd_before^2,
         logratio = 2*(log(sd_after) - log(sd_before))) |> 
  select(contrast, ratio, logratio) |> 
  pivot_longer(everything())

# Mean and 95% credible intervals for contrasts

posterior_2025 |> 
  group_by(name) |> 
  median_qi()

# Posterior probability that sd is lower in exercise 2

posterior_2025 |> 
  mutate(thres = if_else(name == "ratio", 1, 0)) |> 
  group_by(name) |> 
  summarise(mean(value < thres))

#or alternatively
hypothesis(brm_rich_2025, 
           "Detector__Intercept:exercise1^2 > Detector__Intercept:exercise2^2",
           class = "sd")

### plots -------------------------------------------------

#' From this we can see that the prior for SD is not informative, and that the posterior has converged to a tighter probability region.

#prior for sd vs posteriors
brm_rich_2025 |> 
  as_draws_df() |> 
  select(sd_before = 'sd_Detector__Intercept:exercise1',
         sd_after = 'sd_Detector__Intercept:exercise2') |> 
  pivot_longer(1:2, names_to = "sd", values_to = "draw") |> 
  ggplot(aes(x = draw, fill = sd)) +
  stat_slab(aes(xdist = dist_truncated(dist_normal(0, 5), lower = 0)),
            fill = "cornflowerblue", alpha = 0.5, inherit.aes = F,
            data = tibble()) +
  geom_density(color = NA, alpha = 0.5)

#' These are alternative visualization of the contrast between exercises.

#difference

plot_rich_25 <- posterior_2025 |> 
  filter(name == "contrast") |>
  plot_post() +
  labs(title = "2025",
       x = "Standard deviation difference (EX2 - EX1)",
       y = "") 

plot_rich_25

#ratio

posterior_2025 |> 
  filter(name == "ratio") |>
  select(value) |> 
  mcmc_areas(border_size = 1, prob = 0.5, prob_outer = 0.95, point_est = "median")+
  geom_vline(xintercept = 1, linetype = 2)+
  coord_cartesian(xlim = c(0, NA), ylim = c(1,1.1))+
  labs(title = "Posterior detector variance ratio (after / before)",
       x = "Variance ratio", y = "Density")+
  theme_minimal()+
  theme(axis.text.y = element_blank())

#logratio

posterior_2025 |> 
  filter(name == "logratio") |> 
  plot_post()+
  labs(title = "Posterior detector variance log-ratio log(after / before)",
       x = "Variance log-ratio", y = "")


## 2023 ------------------------------------------------------------------

brm_rich_2023 <- brm(formula = rich_formula,
                     data = richness_2023,
                     prior = prior_rich,
                     file = "fits/rich_2023",
                     seed = 967812)

prior_summary(brm_rich_2023)

#posterior summary
summary(brm_rich_2023)

brm_rich_2023 |> 
  mcmc_intervals(pars = vars(b_exercise2:sigma))

posterior_2023 <- as_draws_df(brm_rich_2023) |> 
  select(sd_before = 'sd_Detector__Intercept:exercise1',
         sd_after = 'sd_Detector__Intercept:exercise2') |> 
  mutate(contrast = sd_after - sd_before,
         ratio = sd_after^2/sd_before^2,
         logratio = 2*(log(sd_after) - log(sd_before))) |> 
  select(contrast, ratio, logratio) |> 
  pivot_longer(everything())

# Mean and 95% credible intervals

posterior_2023 |> 
  group_by(name) |> 
  median_qi()

# Posterior probability that sd is lower in exercise 2

posterior_2023 |> 
  mutate(thres = if_else(name == "ratio", 1, 0)) |> 
  group_by(name) |> 
  summarise(mean(value < thres))
  
### plots----------------------------------------------------------------

#prior for sd vs posteriors
brm_rich_2023 |> 
  as_draws_df() |> 
  select(sd_before = 'sd_Detector__Intercept:exercise1',
         sd_after = 'sd_Detector__Intercept:exercise2') |> 
  pivot_longer(1:2, names_to = "sd", values_to = "draw") |> 
  ggplot(aes(x = draw, fill = sd)) +
  stat_slab(aes(xdist = dist_truncated(dist_normal(0, 5), lower = 0)),
            fill = "cornflowerblue", alpha = 0.5, inherit.aes = F,
            data = tibble()) +
  geom_density(color = NA, alpha = 0.5)

#difference

plot_rich_23 <- posterior_2023 |> 
  filter(name == "contrast") |>
  plot_post() +
  labs(title = "2023",
       x = "Standard deviation difference (EX2 - EX1)",
       y = "")

plot_rich_23

#ratio

posterior_2023 |> 
  filter(name == "ratio") |>
  select(value) |> 
  mcmc_areas(border_size = 1, point_est = "median")+
  geom_vline(xintercept = 1, linetype = 2)+
  coord_cartesian(xlim = c(0, 10), ylim = c(1, 1.1))+
  labs(title = "Posterior detector variance ratio (after / before)",
       x = "Variance ratio", y = "Density")+
  theme_minimal()+
  theme(axis.text.y = element_blank())

#logratio

posterior_2023 |> 
  filter(name == "logratio") |>
  plot_post() +
  labs(title = "Detector variance log-ratio log(after / before)",
       x = "Variance log-ratio", y = "")

#' # Models based on dissimilarity
#' Jaccard and Bray-Curtis dissimilarities are bounded between 0 and 1. Euclidean dissimilarity is also bounded in practice, although the upper bound depends on the data. In this case, it is also scaled in the (0, 1) range so that the same models can be used for all three. The likelihood of choice is beta, with logit link to the predictors.

# beta models ------------------------------------------------------------

beta_formula <- brmsformula(dissimilarity ~ exercise + (1|SU) + (1|detector),
                            family = Beta(link_phi = "identity"))

get_prior(beta_formula, data = jaccard_2025)

## priors---------------------------------------------------------------

### visualize beta ----------------------------------------------------
#adapted from https://solomonkurz.netlify.app/blog/2023-06-25-causal-inference-with-beta-regression/

expand_grid(mu = c(.5, .75, .9),
            phi = c(0.1, 1, 10, 100, 500),
            beta = seq(from = .001, to = .999, length.out = 300)) |> 
  mutate(d   = dbeta(beta, shape1 = mu * phi, shape2 = phi - (mu * phi)),
         mu  = str_c("mu==", mu),
         phi = str_c("phi==", phi)) %>% 
  mutate(phi = factor(phi, levels = str_c("phi==", c(0.1, 1, 10, 100, 500)))) |> 
  ggplot(aes(x = beta, y = d)) +
  geom_area(fill = "grey50") +
  scale_x_continuous(breaks = 0:2 / 2, labels = c("0", ".5", "1")) +
  scale_y_continuous("density", breaks = NULL) +
  ggtitle(expression("Beta given different combinations of "*mu*" and "*phi)) +
  coord_cartesian(ylim = c(0, 14)) +
  facet_grid(mu ~ phi, labeller = label_parsed)

### set priors --------------------------------------------------------------

# the inverse logit function
ggplot()+
  geom_function(fun = plogis)+
  xlim(-8, 8)+
  ylab("inv_logit(x)")

#intercept

#' A prior that is too wide results in too much probability mass on extreme values in the outcome scale because of the link function.

# default prior

tibble(x = rstudent_t(1e5, 3, 0, 2.5) |> 
         plogis()) |> 
  ggplot(aes(x = x))+
  stat_slab(fill = "slategray3")+
  ggtitle("Default prior for Intercept on outcome scale")

# adjusted prior

tibble(x = rstudent_t(1e5, 7, 0, 1.5) |> 
         plogis()) |> 
  ggplot(aes(x = x))+
  stat_slab(fill = "rosybrown1")+
  ggtitle("Adjusted prior for Intercept on outcome scale")  

# exercise

tibble(int = rstudent_t(1e5, 7, 0, 1.5),
       b0 = plogis(int),
       b1 = plogis(int + rnorm(1e5, 0, 1)),
       x = abs(b1 - b0)) |> 
  ggplot(aes(x = x))+
  stat_halfeye(fill = "palegreen")+
  ggtitle("Prior for absolute difference of exercise on outcome scale")

#' The default prior for $phi$ puts most of the probability on values smaller than 1, which result in a bimodal distribution. The adjusted prior still allows this, but is concentrated on more reasonable distributions.

# phi

tibble(prior = c("default", "adjusted"),
       dist = c(dist_gamma(0.01, 0.01), dist_gamma(1, 0.02))) |> 
  ggplot(aes(xdist = dist, y = prior)) + 
  stat_halfeye(.width = c(.5, .99)) +
  scale_x_log10(limits = c(0.1, 1000),
                breaks = 10^(-1:3),
                labels = scales::label_number()) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  labs(title = expression("Prior for"~phi))

#check some quantiles

qgamma(0.01, 1, 0.02)
qgamma(0.99, 1, 0.02)

# sd
ggplot()+
  stat_halfeye(aes(xdist = dist_truncated(dist_student_t(3, 0, 2.5), lower = 0)),
               fill = "cornsilk2", .width = c(0.66, 0.95))+
  scale_y_continuous(NULL, breaks = NULL)+
  coord_cartesian(xlim = c(0, 15))+
  ggtitle(expression("Default prior for"~sigma))

#' With this default prior, a SD of 3 is very likely. What does it mean to have a standard deviation of 3 on the outcome scale?

expand_grid(mu = c(-1, 0, 1),
       sd = c(0.5, 1, 3)) |> 
  mutate(x = map2(mu, sd, \(x, y) rnorm(1e5, x, y))) |> 
  unnest(x) |> 
  mutate(mu  = str_c("mu==", mu),
         sd = str_c("sigma==", sd),
         x = plogis(x)) |> 
  ggplot(aes(x = x)) +
  stat_slab(fill = "grey50") +
  scale_y_continuous("density", breaks = NULL) +
  facet_grid(sd ~ mu, labeller = label_parsed) +
  ggtitle("Effect of different random SDs on outcome scale")
  
#' It is clear that SDs larger than 1 should be avoided, since they can shift the intercept to the extremes of the scale, regardless of the mean value. Here is then a more reasonable prior:

ggplot()+
  stat_halfeye(aes(xdist = dist_truncated(dist_student_t(5, 0, 0.5), lower = 0)),
               fill = "cornsilk2", .width = c(0.66, 0.95))+
  scale_y_continuous(NULL, breaks = NULL)+
  ggtitle(expression("Default prior for"~sigma))

prior_beta <- c(prior(student_t(7, 0, 1.5), class = "Intercept"),
                prior(normal(0, 1), class = "b"),
                prior(gamma(1, 0.02), class ="phi", lb = 0),
                prior(student_t(5, 0, 0.5), lb = 0, class = "sd"))

# 2023 has fewer levels, requiring a tighter prior

prior_23 <- c(prior(student_t(7, 0, 1.5), class = "Intercept"),
                prior(normal(0, 1), class = "b"),
                prior(gamma(1, 0.02), class ="phi", lb = 0),
                prior(student_t(3, 0, 0.25), lb = 0, class = "sd"))


data_prior <- data.frame(dissimilarity = 0.5,
                         exercise = factor(1:2),
                         SU = factor(1:2),
                         detector = factor(1:2))

prior_pred <- brm(beta_formula,
                  data = data_prior,
                  sample_prior = "only",
                  prior = prior_beta,
                  file = "fits/beta_prior",
                  seed = 9512)

#' Here's a prior predictive check for beta distributions.

prior_pred |> 
  spread_draws(b_Intercept, b_exercise2, sd_SU__Intercept, sd_detector__Intercept, phi) |> 
  slice_sample(n = 100) |> 
  mutate(su = rnorm(n(), 0, sd_SU__Intercept),
         det = rnorm(n(), 0, sd_detector__Intercept),
         ex1 = plogis(b_Intercept + su + det),
         ex2 = plogis(b_Intercept + b_exercise2, su + det)) |>
  pivot_longer(c(ex1, ex2)) |> 
  expand_grid(x = seq(0, 1, length.out = 100)) |> 
  mutate(density = dbeta(x, value*phi, (1-value)*phi)) |> 
  ggplot(aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.2, color = "navy")+
  facet_wrap(~ name)

## 2025 ------------------------------------------------------

### jaccard -------------------------------------------------------------

brm_jac_2025 <- brm(formula = beta_formula,
                   data = jaccard_2025,
                   prior = prior_beta,
                   file = "fits/jac_2025",
                   seed = 85628)

prior_summary(brm_jac_2025)

summary(brm_jac_2025)

#' It is clear from the following that the prior is very uninformative. This is on the logit scale.

#posterior vs prior for exercise

brm_jac_2025 |>
  post_prior_beta()

#on the logit scale
brm_jac_2025 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Jaccard dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_jac_2025 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Jaccard dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_jac_25 <- brm_jac_2025 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), 
                     limits = c(-0.25, 0.05))+
  labs(title = "Jaccard", x = "Change in dissimilarity (%)", y = "")

plot_lift_jac_25

#quantile dotplot of predictions

plot_pred_jac_25 <- brm_jac_2025 |> 
  plot_pred() +
  labs(x = "Jaccard", y = "Exercise")

plot_pred_jac_25

### bray-curtis ----

brm_bray_2025 <- brm(formula = beta_formula,
                data = bray_2025,
                prior = prior_beta,
                file = "fits/bray_2025",
                seed = 61278)

summary(brm_bray_2025)

#posterior vs. prior

brm_bray_2025 |> 
  post_prior_beta()

#on the logit scale
brm_bray_2025 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Bray-Curtis dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_bray_2025 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Bray-Curtis dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_bray_25 <- brm_bray_2025 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(-.35, .15))+
  labs(title = "Bray-Curtis", x = "Change in dissimilarity (%)", y = "")

plot_lift_bray_25

#quantile dotplot of predictions

plot_pred_bray_25 <- brm_bray_2025 |> 
  plot_pred() +
  labs(x = "Bray-Curtis", y = "Exercise")

plot_pred_bray_25

### euclidean ---------------------------------------------------------

euc_2025 |> 
  ggplot(aes(x = dissimilarity, fill = exercise))+
  geom_histogram(position = "dodge")+
  facet_wrap(~SU)

#model
brm_euc_2025 <- brm(formula = beta_formula,
                    data = euc_2025,
                    prior = prior_beta,
                    file = "fits/euc_2025",
                    seed = 21856)

summary(brm_euc_2025)

prior_summary(brm_euc_2025)

#posterior vs. prior

brm_euc_2025 |> 
  post_prior_beta()

#on the logit scale
brm_euc_2025 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Euclidean dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_euc_2025 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Euclidean dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_euc_25 <- brm_euc_2025 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(-.6, 0))+
  labs(title = "Euclidean", x = "Change in dissimilarity (%)", y = "")

plot_lift_euc_25

#quantile dotplot of predictions

plot_pred_euc_25 <- brm_euc_2025 |> 
  plot_pred() +
  labs(x = "Euclidean", y = "Exercise")

plot_pred_euc_25

## 2023 ---------------------------------------------------

### jaccard ----------------------------------------------------------

brm_jac_2023 <- brm(formula = beta_formula,
                    data = jaccard_2023,
                    prior = prior_23,
                    file = "fits/jac_2023",
                    seed = 59293,
                    control = list(adapt_delta = 0.9))

summary(brm_jac_2023)

prior_summary(brm_jac_2023)

#posterior vs. prior

brm_jac_2023 |> 
  post_prior_beta()

#on the logit scale
brm_jac_2023 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Jaccard dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_jac_2023 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Jaccard dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_jac_23 <- brm_jac_2023 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(-0.25, 0.05))+
  labs(x = "Percent change", y = "")

plot_lift_jac_23

#quantile dotplot of predictions

plot_pred_jac_23 <- brm_jac_2023 |> 
  plot_pred() +
  labs(x = "Jaccard", y = "Exercise")

plot_pred_jac_23

### bray-curtis ------------------------------------------------------

brm_bray_2023 <- brm(formula = beta_formula,
                     data = bray_2023,
                     prior = prior_23,
                     file = "fits/bray_2023",
                     seed = 61278,
                     control = list(adapt_delta = 0.9))

summary(brm_bray_2023)

prior_summary(brm_bray_2023)

#posterior vs. prior

brm_bray_2023 |> 
  post_prior_beta()

#on the logit scale
brm_bray_2023 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Bray-Curtis dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_bray_2023 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Bray-Curtis dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_bray_23 <- brm_bray_2023 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(-.35, .15))+
  labs(x = "Percent change", y = "")

plot_lift_bray_23

#quantile dotplot of predictions

plot_pred_bray_23 <- brm_bray_2023 |> 
  plot_pred() +
  labs(x = "Bray-Curtis", y = "Exercise")

plot_pred_bray_23

### euclidean ---------------------------------------------------------

#model
brm_euc_2023 <- brm(formula = beta_formula,
                    data = euc_2023,
                    prior = prior_23,
                    file = "fits/euc_2023",
                    seed = 21856,
                    control = list(adapt_delta = 0.9))

summary(brm_euc_2023)

prior_summary(brm_euc_2023)

#posterior vs. prior

brm_euc_2023 |> 
  post_prior_beta()

#on the logit scale
brm_euc_2023 |> 
  plot_post(fit = T) +
  labs(title = "Posterior distribution for Euclidean dissimilarity",
       x = "Coefficient of exercise", y = "")

#on the outcome scale

#as difference

brm_euc_2023 |> 
  plot_comp(contrast) +
  labs(title = "Effect of training on Euclidean dissimilarity",
       x = "Difference in dissimilarity", y = "")

#as percent lift
plot_lift_euc_23 <- brm_euc_2023 |> 
  plot_comp(lift) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(-.6, 0))+
  labs(x = "Percent change", y = "")
  
plot_lift_euc_23

#quantile dotplot of predictions

plot_pred_euc_23 <- brm_euc_2023 |> 
  plot_pred() +
  labs(x = "Euclidean", y = "Exercise")

plot_pred_euc_23

# Building the figures -----------------------------------------------------

#richness

wrap_plots(grid::textGrob("Observer richness variability before and after training", 
                          gp = grid::gpar(fontsize = 12)),
           plot_rich_25, plot_rich_23) +
  plot_layout(axes = 'collect',
              ncol = 2, nrow = 2, byrow = T,
              heights = c(0, 1, 1),
              design = "AA
                        BC") +
  plot_annotation(tag_levels = list(c("", paste0(letters[1:2], ")")))) &
  # plot_annotation("Observer richness variability before and after training",
  #                 tag_levels = 'a',
  #                 tag_suffix = ')') 
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        plot.tag = element_text(size = 11),
        plot.tag.location = 'plot',
        margins = margin(t = 10, b = 5, l = 10, r = 5))

ggsave('output/richness.png', bg = 'white', scale = 1.1)

#beta predictions

r1_pred <- wrap_plots(plot_pred_jac_25, plot_pred_bray_25, plot_pred_euc_25,
                      axes = 'collect') &
  theme(axis.title.x = element_blank(),
        margins = margin_part(l = 10, r = 10))

r2_pred <- wrap_plots(plot_pred_jac_23, plot_pred_bray_23, plot_pred_euc_23,
                      axes = 'collect') &
  theme(margins = margin_part(l = 10, r = 10))

wrap_plots(grid::textGrob("Observer dissimilarity in 2025", gp = grid::gpar(fontsize = 12)),
           r1_pred,
           grid::textGrob("Observer dissimilarity in  2023", gp = grid::gpar(fontsize = 12)),
           r2_pred, 
           ncol = 1, heights = rep(c(1, 10), times = 2)) +
  plot_annotation(tag_levels = list(c("", paste0(letters[1:3], ")"),
                                      "", paste0(letters[4:6], ")")))) &
  theme(axis.text = element_text(size = 11),
        plot.tag = element_text(size = 10),
        plot.tag.location = 'panel',
        plot.tag.position = c(-.02, 1))

ggsave('output/beta_pred.png', bg = 'white', scale = 1.1)


#beta lift %

r1_lift <- wrap_plots(plot_lift_jac_25, plot_lift_bray_25, plot_lift_euc_25,
                      axes = 'collect') &
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        margins = margin_part(l = 10, r = 10))

r2_lift <- wrap_plots(plot_lift_jac_23, plot_lift_bray_23, plot_lift_euc_23,
                      axes = 'collect') &
  theme(margins = margin_part(l = 10, r = 10))

wrap_plots(grid::textGrob("Change in dissimilarity in 2025", gp = grid::gpar(fontsize = 12)),
           r1_lift,
           grid::textGrob("Change in dissimilarity in 2023", gp = grid::gpar(fontsize = 12)),
           r2_lift, 
           ncol = 1, heights = rep(c(1, 10), times = 2)) +
  plot_annotation(tag_levels = list(c("", paste0(letters[1:3], ")"),
                                      "", paste0(letters[4:6], ")")))) &
  theme(plot.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        plot.tag = element_text(size = 10),
        plot.tag.location = 'panel',
        plot.tag.position = c(-.02, 1))

ggsave('output/beta_lift.png', bg = 'white', scale = 1.1)

#' # Session info 

sessionInfo()
