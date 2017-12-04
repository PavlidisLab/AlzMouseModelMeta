rm(list=ls())
library(lme4)
str(sleepstudy)
sleepstudy <- sleepstudy
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
summary(fm1)

load("/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model/random_slope_include_NA/early/mixed_model_results.Rdata")
df_slope <- df_all
load("/home/bzhuang/ND_project_combined/AD_mouse_model_project/mixed_model/random_intercept_include_NA/early/mixed_model_results.Rdata")
df_inter <- df_all

i='Fbxo18'
i='Sostdc1'
i='Tmsb10'

df <- as.data.frame(t(array_dat[i, ])) # expression of the genes from all samples

i <- colnames(df)

print(i)
colnames(df) <- 'Expr'
df$Sample <- row.names(df)
df <- noWarnings(left_join(df, array_design))
row.names(df) <- df$Sample
## update level order
df$Disease_stage <- factor(df$Disease_stage, levels = c('WT', 'Disease'))
if('Gender' %in% colnames(df)){ ## set gender M as baseline
    if(length(levels(df$Gender)) ==2){
        df$Gender <- factor(df$Gender, levels = c('M', 'F'))
        gender_col <- 'GenderF'
    }else{gender_col <- NA}
}


## slope
full_model <- lmer(Expr ~ Disease_stage + Gender + (1+ Disease_stage| Study), data = df, REML = F)
null_model = lmer(Expr ~ Gender + (1+Disease_stage|Study), data=df, REML=F)


## intercept
full_model <- lmer(Expr ~ Disease_stage + Gender + (1| Study), data = df, REML = F)
null_model = lmer(Expr ~ Gender + (1|Study), data=df, REML=F)

(anova_r <- anova(full_model,null_model))
(pvalue <- anova_r[, "Pr(>Chisq)"][2])
(chisq <- anova_r[, "Chisq"][2])
(chi_df <- anova_r[, "Chi Df"][2])

fm2 <- lmer(Expr ~ Disease_stage + Gender + (1| Study), data = df, REML = F)

full_model <- fm
(df_t <- data.frame(geneSymbol =i, t(coef(summary(full_model))['Disease_stageDisease',])))


# report all coefficients
## get coefficients for each study
(x <- as.data.frame(t(coefficients(full_model)$Study)))
# get intercepts
(all_intercepts <- as.data.frame(x['(Intercept)', ]))
colnames(all_intercepts) <- paste0(colnames(all_intercepts), '_intercept')
# get fixed effects: (fixed intercept model)
df_effect <- as.data.frame(x["Disease_stageDisease", ])
colnames(df_effect) <- paste0(colnames(df_effect), '_FE_Disease')
if(is.null(gender_col)){ ## if gender is not in the model
    df_effect <-cbind(df_effect, data.frame(FE_Gender = NA, Gender = levels(factor(df$Gender))))
}else{ ## if gender is in the model
    df_effect <-cbind(df_effect, data.frame(FE_Gender = x["GenderF", 1], Gender = gsub('Gender', '', gender_col)))
}


summary(full_model)
x
plot(df$Disease_stage, df$Expr)
plotXY(df,  'Disease_stage', 'Expr', facet_grid = "Study")
