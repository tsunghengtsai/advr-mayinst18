library(tidyverse)
load("adv-R-twin.RData")


# Normalization -----------------------------------------------------------

twin_dia <- tbl_df(twin_dia) %>% 
    mutate(
        log2inty_h = log2(intensity_h), 
        log2inty_l = log2(intensity_l)
    )

med_dia <- twin_dia %>% 
    group_by(run) %>% 
    summarise(log2inty_med = median(log2inty_h, na.rm = TRUE)) %>% 
    mutate(log2inty_adj = median(log2inty_med) - log2inty_med)

twin_dia2 <- left_join(twin_dia, med_dia) %>% 
    mutate(
        log2inty_h = log2inty_h + log2inty_adj, 
        log2inty_l = log2inty_l + log2inty_adj, 
        intensity_h = 2 ^ log2inty_h,
        intensity_l = 2 ^ log2inty_l
    )

twin_dia2


# Linear model to summarize feature log-intensities (one protein) ---------

# One protein (AACT) measured in 3 runs (R001, R002, R003)
sub_dia <- twin_dia2 %>% 
    filter(protein == "AACT", run %in% c("R001", "R002", "R003"))

# Linear model with fixed effects of run and feature, no intercept
fit <- lm(log2inty_l ~ 0 + run + feature, data = sub_dia)

# Extract run effect from the fitted model
summary(fit)
coef(fit)
coef(summary(fit))

fit_coef <- coef(fit)

tibble(effect = names(fit_coef), estimate = fit_coef) %>% 
    filter(grepl("run", effect)) %>% 
    mutate(effect = gsub("run", "", effect))


# Summarization for all proteins with a for loop --------------------------

### Approach 1

### Overview on Rmd

prots <- unique(twin_dia2$protein)  # Proteins in the dataset
df_allprot <- NULL
for (i in seq_along(prots)) {
    oneprot <- twin_dia2 %>% filter(protein == prots[i])  # Subset for one protein
    fit <- lm(log2inty_l ~ 0 + run + feature, data = oneprot)  # Fit a linear model
    fit_coef <- coef(fit)  # Extract fitted parameters
    # Assemble the result in a data frame, and combine with those in previous iterations
    df_oneprot <- tibble(
        protein = prots[i], 
        effect = names(fit_coef), 
        estimate = fit_coef
    ) %>% 
        filter(grepl("run", effect)) %>% 
        mutate(effect = gsub("run", "", effect))
    df_allprot <- rbind(df_allprot, df_oneprot)
}


### Approach 2

### Overview on Rmd

list_allprot <- vector("list", length = length(prots))
for (i in seq_along(prots)) {
    oneprot <- twin_dia2 %>% filter(protein == prots[i])  # Subset for one protein
    fit <- lm(log2inty_l ~ 0 + run + feature, data = oneprot)  # Fit a linear model
    fit_coef <- coef(fit)  # Extract fitted parameters
    # Assemble the result in a data frame, and save to the list
    list_allprot[[i]] <- data_frame(
        protein = prots[i], 
        effect = names(fit_coef), 
        estimate = fit_coef
    ) %>% 
        filter(grepl("run", effect)) %>% 
        mutate(effect = gsub("run", "", effect))
}
df_allprot <- bind_rows(list_allprot)


### Comments on Rmd


# broom -------------------------------------------------------------------

library(broom)
?broom

fit <- lm(log2inty_l ~ 0 + run + feature, data = sub_dia)
summary(fit)
coef(summary(fit))  # Is it tidy?

# Each row is a parameter
tidy(fit)

# Each row is an observation; columns from the fitted model start with .
augment(fit)

# One row for the model
glance(fit)


### Approach 3

### Overview on Rmd

twin_dia2 %>% 
    group_by(protein) %>%
    do(tidy(lm(log2inty_l ~ 0 + run + feature, data = .))) %>% 
    filter(grepl("run", term)) %>% 
    mutate(term = gsub("run", "", term))


# List-columns ------------------------------------------------------------

fit_dia <- twin_dia2 %>% 
    group_by(protein) %>% 
    do(fit = lm(log2inty_l ~ 0 + run + feature, data = .))
fit_dia

# Parameter-level summaries
fit_dia %>% tidy(fit)

# Observation summaries
fit_dia %>% augment(fit)

# Model-level summaries
fit_dia %>% glance(fit)


### Comments on Rmd


# Nested data frame -------------------------------------------------------

nested_dia <- twin_dia2 %>% 
    group_by(protein) %>% 
    nest()
nested_dia

# Same as 
twin_dia2 %>% nest(-protein)

# Data of protein A1AG_BOVINE
nested_dia$data[[1]]

# Same as
nested_dia[[1, "data"]]

nested_dia %>% unnest(data)



### The map functions of purrr

### Overview on Rmd


# For each protein, fit a linear model with purrr::map()
nested_dia <- nested_dia %>% 
    mutate(fit = map(data, ~ lm(log2inty_l ~ 0 + run + feature, data = .)))
nested_dia

# For each fitted model, tidy output object with broom::tidy()
nested_dia <- nested_dia %>% 
    mutate(param = map(fit, tidy))
nested_dia

# Alternative with lapply
nested_dia %>%
    mutate(fit2 = lapply(data, function(x) lm(log2inty_l ~ 0 + run + feature, data = x))) %>% 
    mutate(param2 = lapply(fit2, tidy))


head(nested_dia$param[[1]])

nested_dia %>% 
    mutate(sd_lm = map_dbl(fit, sigma))


# Unnest data back to the original form
nested_dia %>% unnest(param)


### Approach 4

### Overview on Rmd

# Create nested data frame 
nested_dia <- twin_dia2 %>% 
    group_by(protein) %>% 
    nest()

# Model and tidy output objects
nested_dia <- nested_dia %>% 
    mutate(fit = map(data, ~ lm(log2inty_l ~ 0 + run + feature, data = .))) %>% 
    mutate(param = map(fit, tidy))

# Unnest data and additional manipulation
nested_dia %>% unnest(param) %>% 
    filter(grepl("run", term)) %>% 
    mutate(term = gsub("run", "", term))

