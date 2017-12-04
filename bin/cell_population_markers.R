#' 2016-10-07 recorded done the cell type markers for future heatmaps
#' updated 2017-01-16
library('HelperFunctions')
#library("markerGenesManuscript", lib.loc="/space/opt/R-3.1.2/lib64/R/library")
library(markerGeneProfile)

region_ls = c('Hippocampus','Striatum')



## do some quick fixe of the markers 
newcholin = c("Ecel1","Slc18a3","Crabp2","Zic3","Trpc3","Slc5a7","Lhx8") ## 2016-11-01
mouseMarkerGenes$Striatum$ForebrainCholin = newcholin
mouseMarkerGenes$Striatum$Dopaminergic  = mouseMarkerGenes$Midbrain$Dopaminergic    ## some dopaminergic expressed


for (region in region_ls){
    print(region)
    df = eval(parse(text = paste0('mouseMarkerGenes$', region)))
    outdir=paste0('../../doc/cell_type_markers/',Sys.Date(), '/', region, '/')
    dir.create(outdir, showWarnings=F, recursive=T)
    
    for(i in 1:length(df)){
        print(names(df[i]))
        df_out= data.frame( geneSymbols = df[[i]], cell_type = names(df[i]),region = region)
        f_out= paste0(outdir, region, '_', names(df[i]), '.tsv')
        writeTable(df_out, f_out)
    }
    
}

f_out= paste0('../../doc/cell_type_markers/',Sys.Date(), '/mouseMarkerGenes.Rdata')
save(mouseMarkerGenes, file=f_out)
print(f_out)

