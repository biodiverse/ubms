NAME = $(shell grep 'Package:' DESCRIPTION | cut -d ' ' -f2)
VER = $(shell grep 'Version:' DESCRIPTION | cut -d ' ' -f2)

install:
	R CMD INSTALL .

build:
	cd ..; R CMD build $(NAME)

check:
	Rscript -e "devtools::check()"

document:
	Rscript -e "devtools::document()"

README:
	Rscript -e "rmarkdown::render('README.Rmd')"

clean-install:
	R CMD INSTALL --preclean .
