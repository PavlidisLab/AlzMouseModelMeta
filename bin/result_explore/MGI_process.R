# 2016-08-05, process annotations, remove duplicated gene symbols get info from MGI

folder <- '/home/bzhuang/ND_project_combined/ND_results/MGI_info/'
f= paste0(folder, 'annotations.tsv')
fo <- paste0(folder, 'annotation_precessed.tsv')
df <- read.delim(f)
df_tmp <- df[, -c(11,12)]
dup_i <- which(duplicated(df_tmp) == T)
df <- df[-dup_i, ]%>%droplevels()


dup_i <- which(duplicated(df$geneSymbol) == T & df$Input.Type != 'current symbol')
df_t <- df[dup_i, ]%>%droplevels()
df <- df[-dup_i, ]%>%droplevels()


dup_i <- which(duplicated(df$geneSymbol) == T)
dup_genes <- as.character(df$geneSymbol[dup_i])
dup_i <- which(df$geneSymbol %in% dup_genes)

df_part1 <- df[-dup_i, ]%>%droplevels()

df_t <- df[dup_i, ]%>%droplevels()
keep_i <- which(duplicated(df_t$geneSymbol)==T & df_t$Input.Type == 'current symbol')

df_part2 <- df_t[keep_i, ]%>%droplevels()


df_final <- rbind(df_part1, df_part2)
df_final <- orderCol(df_final, cols='geneSymbol')

## change to sense (+) antisense (-)
df_final$Strand <- as.character(df_final$Strand)
index <- which(df_final$Strand == '+')
df_final$Strand[index] <- 'sense'
index <- which(df_final$Strand == '-')
df_final$Strand[index] <- 'antisense'

writeTable(df_final, f_out=fo)
#################################
## get the human disease
f= paste0(folder, 'omim_id.txt')
fo <- paste0(folder, 'annotation_precessed_omim.tsv')
df <- read.delim(f)

df <- df[which(df$OMIM.ID !=""), ] %>% droplevels()

## concatenate:
colnames(df)
df_t <- aggregate(OMIM.ID ~ Input, df, function(x) paste0(x, collapse = '|'))
df_t2 <- aggregate(Term ~ Input, df, function(x) paste0(x, collapse = '|'))
df <- left_join(df_t, df_t2)
colnames(df) <- c('geneSymbol', 'OMIM_ID', 'Human_Disease')


df_omim <- left_join(df_final, df)
writeTable(df_omim, f_out=fo)

#########################
## add Mammalian_Phenotype

f= paste0(folder, 'Mammalian_Phenotype.txt')
fo <- paste0(folder, 'annotation_precessed_omim_pheno.tsv')
df <- read.delim(f)
## with detected
df <- df[which(df$Detected >0), ] %>% droplevels()
df <- df[, c('Input', 'Anatomical.Structure')] %>% droplevels()
colnames(df) <- c('geneSymbol', 'mammalian_phenotype')
df$mammalian_phenotype <- as.character(df$mammalian_phenotype)
df$mammalian_phenotype2 <- unlist(strsplit(df$mammalian_phenotype, split = ': '))

df$mammalian_phenotype2 <- sapply(df$mammalian_phenotype,  function(x) unlist(strsplit(x, split = ": "))[2])

df2 <- aggregate(mammalian_phenotype2 ~ geneSymbol, df, function(x) paste0(x, collapse = '|'))

colnames(df2)[2] <- 'mammalian_phenotype'


df_omim_phenotype <- left_join(df_omim, df2)
writeTable(df_omim_phenotype, f_out=fo)
