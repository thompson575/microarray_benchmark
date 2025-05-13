library(Biobase)
library(GEOquery)
library(fs)

# Home directory
home = "C:/Projects/RDevelopment/MLforGE/microarray_benchmark"

# Read in gse format
# gse <- getGEO( file = path(home, "data/rawData/GSE252168_series_matrix.txt.gz"))

# Save in the temp folder
# saveRDS(gse, path(home, "temp/GSE252168.rds"))

# Read from the temp folder
#gse <- readRDS(path(home, "temp/GSE252168.rds"))

#features <- as.data.frame(pData(gse))

#saveRDS(features, path(home, "temp/features_df.rds") )
# Transpose the expression measures
# df <- as.data.frame( t(exprs(gse)))
# saveRDS(df, path(home, "temp/GSE252168_df.rds"))

# df has 303 rows (patients) and 30715 columns (genes)
#  180 controls and 123 Lung Cancer

# Read df from temp
df <- readRDS(path(home, "temp/GSE252168_df.rds"))
df <- log2(df)
features <- readRDS(path(home, "temp/features_df.rds") )
df$y <- factor(features[["status:ch1"]])
df$sex <- factor(features[["gender:ch1"]])
df$biobank <- factor(features[["biobank:ch1"]])
# x <- log2(df[, 5])
# par(mfrow=c(2, 1))
# hist(x[df$y == "case"], breaks=50, main="Lung Ca", xlab="Expression",
#      xlim=c(quantile(x, 0.05), quantile(x, 0.95)))
# hist(x[df$y == "control"], breaks=50, main="control", xlab="Expression", 
#      xlim=c(quantile(x, 0.05), quantile(x, 0.95)))
# par(mfrow=c(1, 1))

x <- df[, 13764]
boxplot(x ~ df$y)

summary(x[df$y == "case"])
summary(x[df$y == "control"])

f <- function(x) c(m=mean(x), sd=sd(x))
aggregate( x, by = list(df$y), FUN=f)

p <- rep(0, 30715)
for(i in 1:30715) {
  p[i] <- t.test(df[, i] ~ df$y)$p.value
}
which(p < 1e-50)
min(p)
p <- -log10(sort(p))
q <- -log10( (1:30715-0.5)/30715)
plot(q, p, pch=16, xlab="Expected -log10(p)", ylab="-log10(p)")
abline(a=0, b=1)

# 13764 14520 21923 28849
plot(df[,13764], df[,21923], type="n")
points(df[df$y=="case" ,13764], df[df$y=="case",21823], col="red")
points(df[df$y=="control" ,13764], df[df$y=="control",21923], col="blue")

pairs(df[, c(13764,14520,21923,28849)], col=c("red","blue")[df$y])

table(df$y, df$sex)
p <- rep(0, 30715)
for(i in 1:30715) {
  p[i] <- t.test(df[, i] ~ df$sex)$p.value
}
p <- -log10(sort(p))
q <- -log10( (1:30715-0.5)/30715)
plot(q, p, pch=16, xlab="Expected -log10(p)", ylab="-log10(p)")
abline(a=0, b=1)

pairs(df[, c(13764,14520,21923,28849)], col=c("red","blue")[df$sex])

pairs(df[, c(13764,14520,21923,28849)], col=c("green","red","blue")[df$biobank])

# ===============================================================
# Adjustment for Blood tube and Gender
#
df$tube <- factor( df$biobank == "HUNT", labels=c("PAXgene", "tempus"))
for(i in 1:30715) {
  mod <- lm(df[, i] ~ df$tube + df$sex)
  df[, i] <- mod$residuals
}
saveRDS(df, path(home, "temp/GSE252168_Adjusted_df.rds"))

df <- readRDS(path(home, "temp/GSE252168_Adjusted_df.rds"))
p <- rep(0, 30715)
for(i in 1:30715) {
  p[i] <- t.test(df[, i] ~ df$y)$p.value
}

p <- -log10(p)
which(p>13.5)
q <- -log10( (1:30715-0.5)/30715)
plot(q, sort(p), pch=16, xlab="Expected -log10(p)", ylab="-log10(p)")
abline(a=0, b=1)

pairs(df[, c(13764,14520,21923,28849)], col=c("red","blue")[df$y])

pairs(df[, c(3044,13328,13764,21923)], col=c("red","blue")[df$y])

# ==================================================================
# mlr3
#
library(mlr3verse)
library(fs)
library(mlr3tuning)

# Home directory
home = "C:/Projects/RDevelopment/MLforGE/microarray_benchmark"

df <- readRDS(path(home, "temp/GSE252168_Adjusted_df.rds"))
df$sex = df$tube = df$biobank = NULL

# ---------------------------------------
# Baseline model
#
Tk=TaskClassif$new(id="GSE252168", df, target="y")
Lr=lrn("classif.simple_logistic")
Ms=msr("classif.acc")
Fi=flt("kruskal_test")

Gr = po("filter", Fi, filter.nfeat = 10) %>>%
  po("learner", Lr)

Lg = as_learner(Gr)

Lg$train(Tk)
# accuracy = 79.5%
Lg$predict(Tk)$score(Ms)

cv = rsmp("cv", folds=5)
cv$instantiate(Tk)

rr = resample(Tk, Lg, cv)
# accuracy = 73.6% (nfeat=30) 72.3%  (nfeat=20)  72.9 (nfeat=10)
rr$aggregate(Ms)

Gr = po("filter", Fi, filter.nfeat = to_tune(5, 50)) %>>%
  po("learner", Lr)

Lg = as_learner(Gr)

instance = ti(
  task = Tk,
  learner = Lg,
  resampling = rsmp("cv", folds = 5),
  measures = Ms,
  terminator = trm("none")
)

Tn = tnr("grid_search", resolution = 7)
Tn$optimize(instance)

as.data.table(instance$archive)[1:7, .(kruskal_test.filter.nfeat, classif.acc)]

library(tidyverse)
as.data.table(instance$archive)[1:7, .(kruskal_test.filter.nfeat, classif.acc)] |>
  as_tibble() |>
  ggplot( aes(x=kruskal_test.filter.nfeat, y=classif.acc)) + geom_point() +
  geom_line()

Fi$calculate(Tk)
Order = as.data.table(Fi) |> as_tibble()
Order |> arrange(desc(score) ) 

hist(Order$score[1:200], breaks=50)
pairs(df[, c("ILMN_1762769","ILMN_3200018",
             "ILMN_1898691","ILMN_1758728")], col=c("red","blue")[df$y])

topRNAs = Order$feature[1:250]
df <- df[, c(topRNAs, "y")]

saveRDS(df, path(home, "temp/GSE252168_Adjusted_top250_df.rds"))

# ==========================================================
# Analyse the top 250 RNAs
#
df <- readRDS(path(home, "temp/GSE252168_Adjusted_top250_df.rds"))
Tk=TaskClassif$new(id="GSE252168", df, target="y")

# -------------------------------------------------
# Logistic
#
Lr=lrn("classif.simple_logistic")
Ms=msr("classif.acc")
Fi=flt("kruskal_test")

Gr = po("filter", Fi, filter.nfeat = 5) %>>%
  po("learner", Lr)
Lg = as_learner(Gr)
cv = rsmp("cv", folds=10)
cv$instantiate(Tk)

rr = resample(Tk, Lg, cv)
rr$aggregate(Ms)

# -------------------------------------------------
# Tree
#
Lr=lrn("classif.rpart")
Lr$param_set$values=list(cp=0.01)
Ms=msr("classif.acc")
cv = rsmp("cv", folds=10)
cv$instantiate(Tk)
rr = resample(Tk, Lr, cv)
rr$aggregate(Ms)

# -------------------------------------------------
# XGBoost
#
Lr=lrn("classif.xgboost")
Lr$param_set$values=list(eta=0.3, nrounds=100)
Ms=msr("classif.acc")
cv = rsmp("cv", folds=10)
cv$instantiate(Tk)
rr = resample(Tk, Lr, cv)
rr$aggregate(Ms)

# -------------------------------------------------
# ranger
#
Lr=lrn("classif.ranger")
Lr$param_set$values=list(mtry=50, num.trees=1000)
Ms=msr("classif.acc")
cv = rsmp("cv", folds=10)
cv$instantiate(Tk)
rr = resample(Tk, Lr, cv)
rr$aggregate(Ms)