# script runs R scripts to: 1) estimate annual averages, 2) make UK-pls predictions, 3) evaluate models

# cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a

## run rscripts
rscript 1_act_annual.R
rscript 2_act_uk.R
rscript 3_model_eval.R

## knit markdown with results 
Rscript -e 'rmarkdown::render("4_results_summary.Rmd", "html_document")'

echo "done running scripts"
