load("adv-R-twin.RData")

str(twin_dia)

head(twin_srm)

library(tidyverse)


# Helpful conventions for data wrangling ----------------------------------

# Creating "tibble"
twin_dia <- as_tibble(twin_dia)
twin_srm <- as_tibble(twin_srm)

# Tibbles are special data frames
class(twin_dia)

# With better print method
twin_dia

# Call RStudio viewer
View(twin_dia)

### The pipe operator
# The pipe operator %>% to chain multiple operations
# The RStudio keyboard shortcut: Ctrl + Shift + M (Windows), Cmd + Shift + M (Mac).


# Tidy data ---------------------------------------------------------------

### Overview on Rmd


# load("messy.RData")

# Making untidy data
td_sub <- as_tibble(twin_srm) %>% select(protein, feature, run, intensity_h, intensity_l) %>% 
    filter(run %in% c("R001", "R002", "R003"), protein %in% c("APOA", "C1QA"))

sub1a <- td_sub %>% select(-intensity_l) %>% spread(run, intensity_h, convert = F)
sub1b <- td_sub %>% select(-intensity_h) %>% spread(run, intensity_l, convert = F)
sub2 <- td_sub %>% gather(key = "label", value = "intensity", intensity_h, intensity_l) %>% 
    mutate(label = ifelse(label == "intensity_h", "heavy", "light")) %>% 
    arrange(protein, feature, run, label)
sub12 <- sub2 %>% spread(key = run, value = intensity)
sub3 <- td_sub %>% unite(intensity_both, intensity_h, intensity_l, sep = "/")


# Some of the column names (R001, R002, R003) are values of a variable, rather than variables
# Subset for heavy channel
sub1a
# Subset for light channel
sub1b

# Different views on sub2
sub2

# A mix of both issues
sub12

# Two variables in one column; values saved as strings
sub3


# Using tidyr -------------------------------------------------------------

### Overview on Rmd

### Use gather() to gather multiple columns into a key-value pair
sub1a

gather(sub1a, key = run, value = intensity_h, R001, R002, R003)
sub1a %>% gather(key = run, value = intensity_h, R001, R002, R003)

# Apply for both sub1a and sub1b and merge the results
tidy1a <- sub1a %>% gather(key = run, value = intensity_h, R001, R002, R003)
tidy1b <- sub1b %>% gather(key = run, value = intensity_l, R001, R002, R003)
dplyr::left_join(tidy1a, tidy1b)  # merge two parts of the dataset, introduced later 


### Use spread() to spread a key-value pair into multiple columns
sub2

spread(sub2, key = label, value = intensity)
sub2 %>% spread(key = label, value = intensity)


### Sometimes you need both gather() and spread() [make it a challenge?]
sub12

sub12 %>% gather(run, intensity, R001, R002, R003)
sub12 %>% 
    gather(run, intensity, R001, R002, R003) %>% 
    spread(key = label, value = intensity) %>% 
    dplyr::rename(intensity_h = heavy, intensity_l = light)


## Use separate() to split a column by a string separator
sub3 %>% separate(col = intensity_both, into = c("intensity_h", "intensity_l"), sep = "/")

# Try to convert to better types using convert = TRUE
sub3 %>% separate(col = intensity_both, into = c("intensity_h", "intensity_l"), sep = "/", convert = TRUE)

# Separate intensity_both and feature columns
sub3 %>% 
    separate(col = intensity_both, into = c("intensity_h", "intensity_l"), sep = "/", convert = TRUE) %>%
    separate(col = feature, into = c("peptide", "z1", "fragment", "z3"), sep = "_")


### Use unite() to merge columns into a single column
sub3 %>% 
    separate(col = intensity_both, into = c("intensity_h", "intensity_l"), sep = "/", convert = TRUE) %>%
    separate(col = feature, into = c("peptide", "z1", "fragment", "z3"), sep = "_") %>% 
    unite(col = transition, z1, fragment, z3, sep = "_")


# Using dplyr -------------------------------------------------------------

### Overview on Rmd


### Use rename() to rename columns
twin_dia %>% rename(inty_H = intensity_h, inty_L = intensity_l)


### Use arrange() to order rows
# Order rows by values of columns protein, run, and feature
twin_dia %>% arrange(protein, run, feature)

twin_dia %>% arrange(desc(subject))


### Use select() to extract existing variables
# Select columns protein and feature
twin_dia %>% select(protein, feature)

# Exclude column pair
twin_dia %>% select(-pair)

# Select from column feature to column intensity_h
twin_dia %>% select(feature:intensity_h)

# Helpful to obtain unique values for particular variables
twin_dia %>% 
    select(protein, feature) %>% 
    distinct()

# Same as
twin_dia %>% distinct(protein, feature)


### Use filter() to extract existing observations
twin_dia %>% filter(!is.na(intensity_h))

# Comma as AND operation
twin_dia %>% filter(is.na(intensity_h), !is.na(intensity_l))


### Use mutate() to add new variables with window functions (vector -> vector)
# Log2 transformation
twin_dia %>% mutate(log2inty_l = log2(intensity_l))

# Use the just generated variables
twin_dia %>% 
    mutate(
        log2inty_h = log2(intensity_h), 
        log2inty_l = log2(intensity_l), 
        log2inty_d = log2inty_l - log2inty_h
    )


### Use group_by() and summarise() to make grouped summaries
###  - summarise() uses summary functions (vector -> single value)
###  - group_by() defines the unit of analysis
# Compute mean, sd and median of values in column intensity_l
twin_dia %>% 
    summarise(
        intensity_ave = mean(intensity_l, na.rm = TRUE), 
        intensity_sd = sd(intensity_l, na.rm = TRUE), 
        intensity_med = median(intensity_l, na.rm = TRUE)
    )

# Compute mean, sd and median of values in column intensity_l, within each run
twin_dia %>% 
    group_by(run) %>% 
    summarise(
        intensity_ave = mean(intensity_l, na.rm = TRUE), 
        intensity_sd = sd(intensity_l, na.rm = TRUE), 
        intensity_med = median(intensity_l, na.rm = TRUE)
    )


### Exercise 2: compute the quantities for constant normalization



# To address Task 1, we will then need to merge this summary back to the original data frame.


# Merge datasets ----------------------------------------------------------

### Overview on Rmd

x <- tibble(
    key = c(1, 2, 3), 
    val_x = c("x1", "x2", "x3")
)

y <- tibble(
    key = c(1, 2, 4), 
    val_y = c("y1", "y2", "y3")
)


### Mutating joins
# Inner join
inner_join(x, y)
x %>% inner_join(y)

# Outer joins
left_join(x, y)
x %>% left_join(y)

right_join(x, y)
x %>% right_join(y)

full_join(x, y)
x %>% full_join(y)

# merge() in base R
merge(x, y)
merge(x, y, all.x = TRUE)
merge(x, y, all.y = TRUE)
merge(x, y, all.x = TRUE, all.y = TRUE)


### Filtering joins
# keep all observations in x that have a match in y
semi_join(x, y)

# Drops all observations in x that have a match in y
anti_join(x, y)


# Task 1: constant normalization ------------------------------------------
# Use summarise() and group_by() to compute the run-level adjustment
twin_dia <- twin_dia %>% 
    mutate(
        log2inty_h = log2(intensity_h), 
        log2inty_l = log2(intensity_l)
    )

med_dia <- twin_dia %>% 
    group_by(run) %>% 
    summarise(log2inty_med = median(log2inty_h, na.rm = TRUE)) %>% 
    mutate(log2inty_adj = median(log2inty_med) - log2inty_med)

med_dia

# Merge the adjustment back to the original dataset
left_join(twin_dia, med_dia)

twin_dia2 <- left_join(twin_dia, med_dia) %>% 
    mutate(
        log2inty_h = log2inty_h + log2inty_adj, 
        log2inty_l = log2inty_l + log2inty_adj, 
        intensity_h = 2 ^ log2inty_h,
        intensity_l = 2 ^ log2inty_l
    )

head(tapply(twin_dia2$log2inty_h, twin_dia2$run, median, na.rm = TRUE))

# Similarly, for the SRM dataset
twin_srm <- twin_srm %>% 
    mutate(
        log2inty_h = log2(intensity_h), 
        log2inty_l = log2(intensity_l)
    ) 

med_srm <- twin_srm %>% group_by(run) %>% 
    summarise(log2inty_med = median(log2inty_h, na.rm = TRUE)) %>% 
    mutate(log2inty_adj = median(log2inty_med) - log2inty_med)

twin_srm2 <- left_join(twin_srm, med_srm) %>% 
    mutate(
        log2inty_h = log2inty_h + log2inty_adj, 
        log2inty_l = log2inty_l + log2inty_adj, 
        intensity_h = 2 ^ log2inty_h,
        intensity_l = 2 ^ log2inty_l
    )

head(tapply(twin_srm2$log2inty_h, twin_srm2$run, median, na.rm=T))


### Visualize the result
# Boxplot of feature log-intensities in each run, before normalization
twin_dia %>% filter(grepl(paste(sprintf("%03d", 1:20), collapse = "|"), run)) %>% 
    ggplot(aes(run, log2inty_h)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Boxplot of feature log-intensities in each run, after normalization
twin_dia2 %>% filter(grepl(paste(sprintf("%03d", 1:20), collapse = "|"), run)) %>% 
    ggplot(aes(run, log2inty_h)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Task 2: comparison between DIA and SRM datasets -------------------------
# Summarization in each dataset
los_dia <- twin_dia2 %>% 
    group_by(run, protein) %>% 
    summarise(sum_dia = sum(intensity_l, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(logsum_dia = ifelse(sum_dia == 0, 0, log2(sum_dia)))
los_dia

# Summarization for the SRM data
los_srm <- twin_srm2 %>% 
    group_by(run, protein) %>% 
    summarise(sum_srm = sum(intensity_l, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(logsum_srm = ifelse(sum_srm == 0, 0, log2(sum_srm)))

# Merge results (with proteins quantified in both)
los_all <- inner_join(los_dia, los_srm)

ggplot(los_all, aes(logsum_dia, logsum_srm)) + 
    geom_point() + geom_smooth(se = FALSE, method = "lm")

ggplot(los_all, aes(logsum_dia, logsum_srm)) + 
    geom_point(aes(colour = protein))


# Compute the correlation coefficient
cor(los_all$logsum_dia, los_all$logsum_srm)

# Compute the correlation per protein
los_all %>% group_by(protein) %>% 
    summarise(correlation = cor(logsum_dia, logsum_srm))


# Apply arbitrary operations to grouped data ------------------------------
cor.test(los_all$logsum_dia, los_all$logsum_srm)

# This would fail... 
los_all %>% group_by(protein) %>% 
    summarise(corres = cor.test(logsum_dia, logsum_srm))

### group_by() + do()
# Here the operation head() returns multiple rows
twin_dia2 %>% 
    group_by(protein) %>% 
    do(head(., 2))

# The pronoun . is used as an argument placeholder, referring to the group data to be processed

# If you use a named argument inside do(), it creates a list-column in the output
twin_dia2 %>% 
    group_by(protein) %>% 
    do(top2 = head(., 2))

# The list-column is useful to store arbitrary R objects, such as models
los_all %>% group_by(protein) %>% 
    do(fit_cor = cor.test(.$logsum_dia, .$logsum_srm))

# We can use double brackets [[]] to retrieve the model objects from the list-column fit_cor
los_cor <- los_all %>% group_by(protein) %>% 
    do(fit_cor = cor.test(.$logsum_dia, .$logsum_srm))
los_cor$fit_cor[[1]]


### More on list-columns in the next section

