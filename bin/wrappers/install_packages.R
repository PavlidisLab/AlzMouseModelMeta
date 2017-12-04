# install required packages

## CRAN packages:
cran_pack <- c('plyr', 'dplyr','ggplot2', 'RColorBrewer', 'reshape2', 'influence.ME','dendextend','tidyr','MetaDE', 'lme4')

(cran_pack=setdiff(cran_pack, rownames(installed.packages())))

for(i in cran_pack){
    install.packages(i, dependencies = T)
}


## biocLite packages:
bio_packs=c('impute','genefilter','limma','globaltest','preprocessCore','simpleaffy','affy','affyPLM','GEOquery',
            'globaltest', 'oligo', 'lumi', "pd.moex.1.0.st.v1")


(bio_packs=setdiff(bio_packs, rownames(installed.packages())))
source('http://bioconductor.org/biocLite.R')
for(i in bio_packs){
    biocLite(i)
}


## install
## see https://github.com/oganm/markerGeneProfile

devtools::install_github('oganm/markerGeneProfile')
