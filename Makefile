NAME = $(shell grep 'Package:' DESCRIPTION | cut -d ' ' -f2)
VER = $(shell grep 'Version:' DESCRIPTION | cut -d ' ' -f2)

install:
	R CMD INSTALL .

build:
	Rscript -e "Rcpp::compileAttributes()"
	cd ..; R CMD build $(NAME)

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

document:
	Rscript -e "devtools::document()"

README:
	Rscript -e "rmarkdown::render('README.Rmd')"
	pandoc README.md -o README.html
	firefox README.html
	sleep 3
	rm README.html

clean-install:
	R CMD INSTALL --preclean .
