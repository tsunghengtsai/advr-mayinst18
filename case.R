load("adv-R-twin.RData")

ls()

str(twin_dia)
str(twin_srm)

head(twin_dia)

class(twin_dia)

# Dimension of the data
dim(twin_dia)
nrow(twin_dia)
ncol(twin_dia)

# Level of the categorical variable
levels(twin_dia$zygosity)

# Call the data viewer in RStudio.
View(twin_dia)


### Data analysis tasks for the case study


# Task 1: median of log-intensities (heavy channel) per run ---------------


### Median of log-intensities in a run 


### Data manipulation and transformation
# Column names
colnames(twin_dia)

# Extract columns for sample annotation
design <- unique(twin_dia[, c("run", "pair", "zygosity", "subject", "visit")])
head(design)

# Filter based on columns intensity_h and run
sub_dia <- twin_dia[!is.na(twin_dia$intensity_h), ]
head(sub_dia)

# Use logical operations to add addition criteria
sub_dia <- twin_dia[!is.na(twin_dia$intensity_h) & twin_dia$run == "R001", ]
head(sub_dia)


# Add new variables with log2 transformation
twin_dia$log2inty_h <- log2(twin_dia$intensity_h)
twin_dia$log2inty_l <- log2(twin_dia$intensity_l)

head(twin_dia)

## Combine multiple operatons 
# Compute the median of feature intensities from the heavy channel (`log2inty_h`) in one run (`R001`):
median(twin_dia$log2inty_h[!is.na(twin_dia$log2inty_h) & twin_dia$run == "R001"])

### Repeating name of data frames 
# with(twin_dia, median(log2inty_h[!is.na(log2inty_h) & run == "R001"]))

### Nested representation


### Median of log-intensities in every run (grouped summaries)


### Approach 1a - for loop
runs <- unique(twin_dia$run)  # unique run ids in the dataset
medians <- rep(0, length(runs))  # create a vector to restore the result
for (i in seq_along(runs)) {
    medians[i] <- median(twin_dia$log2inty_h[twin_dia$run == runs[i]], na.rm = TRUE)
}
str(medians)


### Approach 1b - tapply()
medians <- tapply(twin_dia$log2inty_h, twin_dia$run, median, na.rm = TRUE)
head(medians)  # a named vector is returned


# Approach 1c - aggregate()
df_median <- aggregate(twin_dia$log2inty_h, list(run = twin_dia$run), median, na.rm = TRUE)
head(df_median)  # a data frame is returned

# Rename second column for later use
colnames(df_median)[2] <- "run_median"
head(df_median)


# Task 2: normalization ---------------------------------------------------
# Adjust the log-intensities to equalize the medians across runs to a global median

(gbl_median <- median(medians, na.rm = TRUE))

head(twin_dia)


### Approach 2a: use a for loop
for (ii in names(medians)) {
    log2inty_h <- twin_dia$log2inty_h[twin_dia$run == ii]
    log2inty_l <- twin_dia$log2inty_l[twin_dia$run == ii]
    twin_dia$log2inty_h[twin_dia$run == ii] <- log2inty_h - medians[ii] + gbl_median
    twin_dia$log2inty_l[twin_dia$run == ii] <- log2inty_l - medians[ii] + gbl_median
}

# Check the normalized medians
head(tapply(twin_dia$log2inty_h, twin_dia$run, median, na.rm = TRUE))

# Back to unnormalized values for Approach 2b
twin_dia$log2inty_h <- log2(twin_dia$intensity_h)
twin_dia$log2inty_l <- log2(twin_dia$intensity_l)


### Approach 2b: use a vectorized representation 

# This gives the medians of runs R001, R002, R003, R001
medians[c("R001", "R002", "R003", "R001")]

# Use vectorized representation to normalize the dataset
twin_dia$log2inty_h <- twin_dia$log2inty_h - medians[twin_dia$run] + gbl_median
twin_dia$log2inty_l <- twin_dia$log2inty_l - medians[twin_dia$run] + gbl_median

# Check the normalized medians
head(tapply(twin_dia$log2inty_h, twin_dia$run, median, na.rm = TRUE))

# Back to unnormalized values for Approach 2c
twin_dia$log2inty_h <- log2(twin_dia$intensity_h)
twin_dia$log2inty_l <- log2(twin_dia$intensity_l)


### Approach 2c: merge computed medians to original data frame

twin_dia2 <- merge(x = twin_dia, y = df_median)
head(twin_dia2)

twin_dia2$log2inty_h <- twin_dia2$log2inty_h - twin_dia2$run_median + gbl_median
twin_dia2$log2inty_l <- twin_dia2$log2inty_l - twin_dia2$run_median + gbl_median

head(tapply(twin_dia2$log2inty_h, twin_dia2$run, median, na.rm = TRUE))


# Task 3: summarization of feature intensities ----------------------------

### Log of sum for each protein in each run 


# Transform the normalized log-intensities back to the original scale
twin_dia2$intensity_h <- 2 ^ (twin_dia2$log2inty_h)
twin_dia2$intensity_l <- 2 ^ (twin_dia2$log2inty_l)


### Approach 3a: use two for loops (skipped here)


### Approach 3b: use tapply() with grouping variables defined in a list
sum_t <- tapply(
    twin_dia2$intensity_l, 
    list(run = twin_dia2$run, protein = twin_dia2$protein), 
    function(x) log2(sum(x, na.rm = TRUE))
)
head(sum_t)


### Approach 3c: use aggregate()
sum_g <- aggregate(
    twin_dia2$intensity_l,
    list(run = twin_dia2$run, protein = twin_dia2$protein),
    function(x) log2(sum(x, na.rm = TRUE))
)
head(sum_g)

# Alternatively, use a formula representation
sum_g <- aggregate(
    intensity_l ~ run + protein, 
    data = twin_dia2, 
    function(x) log2(sum(x, na.rm = TRUE))
)
head(sum_g)

# Rename the third column
colnames(sum_g)[3] <- "log2inty"


### matrix vs. data frame for subsequent analysis 


# Task 4: fit a linear model to characterize summarized intensities -------

# For each protein, fit a linear model, extract summaries of the fitted model, 
# and draw model-based inference


## Fit a linear regression model with lm()

# Merge the design information
head(design)
df_sum <- merge(sum_g, design)

head(df_sum)

sub_sum <- df_sum[df_sum$protein == "A2MG", ]  # Subset for protein A2MG
fit <- lm(log2inty ~ zygosity, data = sub_sum)

# Same as
fit <- lm(log2inty ~ zygosity, data = df_sum, subset = df_sum$protein == "A2MG")

class(fit)


## Utility functions

# Use summary() to display an overview of the fitted model
summary(fit)


# Use coef() to retrieve estimated coefficients
coef(fit)
coef(summary(fit))


### Use fitted() to retrieve fitted values
head(fitted(fit))


# More detail
# str(summary(fit))


# To derive summaries for all the proteins, additional efforts are required to 
# extract relevant information and package them into a convenient format. If we 
# want to avoid using for loops, we will need to define specialized functions to 
# extract summaries of interest that can be passed on to tapply() or aggregate()


# Task 5: model-based inference -------------------------------------------

# Two-sample $t$-test on the subset of summaries for protein A2MG
ttest <- t.test(log2inty ~ zygosity, data = sub_sum)
ttest

# Same as 
t.test(log2inty ~ zygosity, data = df_sum, subset = df_sum$protein == "A2MG")

# Alternatively, you can use two vectors for the two samples
t.test(x = sub_sum$log2inty[sub_sum$zygosity == "DZ"], y = sub_sum$log2inty[sub_sum$zygosity == "MZ"])

str(ttest)

ttest$estimate
ttest$statistic
ttest$p.value

# As in Task 4, it requires additional efforts to work with model objects, in 
# order to make and combine summaries for all the proteins. How would you 
# implement a workflow to summarize the means of two groups, difference, 
# t-statistic, p-value in every protein? Where does the inconvenience come from?

