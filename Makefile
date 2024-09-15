.PHONY: vignettes
NAME = $(shell grep 'Package:' DESCRIPTION | cut -d ' ' -f2)
VER = $(shell grep 'Version:' DESCRIPTION | cut -d ' ' -f2)

install:
	R CMD INSTALL .

build:
	Rscript -e "Rcpp::compileAttributes()"
	make document
	cd ..; R CMD build $(NAME)

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "options(testthat.progress.max_fails=30); devtools::test()"

document:
	Rscript -e "devtools::document()"

vignettes: vignettes/random-effects.Rmd vignettes/JAGS-comparison.Rmd vignettes/spatial-models.Rmd vignettes/grouse-example.Rmd
	touch test.png; mv *.png vignettes; rm vignettes/test.png
	Rscript -e "devtools::build_vignettes()"

site:
	Rscript -e "pkgdown::build_site()"
	firefox docs/index.html

clean:
	rm -r src/*.o
	rm -r src/*.so

clean-install:
	R CMD INSTALL --preclean .

coverage:
	Rscript -e \
		'Sys.setenv(NOT_CRAN="true"); covr::report(file="/tmp/ubms-report.html")'
	firefox /tmp/ubms-report.html

vignettes/random-effects.Rmd: vignettes/random-effects.Rmd.orig
	Rscript -e "knitr::knit('vignettes/random-effects.Rmd.orig', output='vignettes/random-effects.Rmd')"

vignettes/JAGS-comparison.Rmd: vignettes/JAGS-comparison.Rmd.orig
	Rscript -e "knitr::knit('vignettes/JAGS-comparison.Rmd.orig', output='vignettes/JAGS-comparison.Rmd')"

vignettes/spatial-models.Rmd: vignettes/spatial-models.Rmd.orig
	Rscript -e "knitr::knit('vignettes/spatial-models.Rmd.orig', output='vignettes/spatial-models.Rmd')"

vignettes/grouse-example.Rmd: vignettes/grouse-example.Rmd.orig
	Rscript -e "knitr::knit('vignettes/grouse-example.Rmd.orig', output='vignettes/grouse-example.Rmd')"
