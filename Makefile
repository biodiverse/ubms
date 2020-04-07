document:
	Rscript -e "devtools::document()"

README:
	Rscript -e "rmarkdown::render('README.Rmd')"
