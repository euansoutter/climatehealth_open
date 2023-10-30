# Run lines below to run API server
devtools::load_all()
library(plumber)

plumber::pr("climatehealth_api.R") %>%
  plumber::pr_run(port=8000)
