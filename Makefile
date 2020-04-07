document:
	Rscript -e "devtools::document()"

README:
	Rscript -e "rmarkdown::render('README.Rmd')"

install:
	R CMD INSTALL .

clean-install:
	R CMD INSTALL --preclean .
