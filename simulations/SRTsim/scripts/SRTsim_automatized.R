library(rmarkdown)
library('SRTsim')

sample_name = 'Visium_FFPE_V43T08-041_D'
preprocessing_output_dir <- './filtered_data/'
sim_output_dir <- '../data/'

#count thresholds for filtering
lower_threshold = 0
upper_threshold = 50000

#threshold to define spots
basal_Ucell_threshold = 0.25
classical_Ucell_threshold = 0.25

#simulation seed (SRTsim_same_domains)
seeds <- c(1,2,3)

print(paste0('#### ',sample_name, ' - Simulation ####'))

render('../../../analysis/spatial_preprocessing.Rmd', 
       output_dir = "../../../analysis/report/", 
       output_file = paste0("report_spatial_filtering_", sample_name,'.html')
)

simSRT  <- createSRT(count_in=filtered_count,loc_in=loc)

for(seed in seeds){
    print(paste(seed,'/',length(seeds)))
    render('SRTsim_same_domains.Rmd',
           output_dir = "../report/",
           output_file = paste0("report_SRTsim_same_domains_", sample_name, '.html')
    )
}

