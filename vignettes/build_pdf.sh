#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pbdDMAT-guide.Rnw
bibtex pbdDMAT-guide
#makeindex pbdDMAT-guide
#pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
Rscript -e "tools::compactPDF('.', gs_quality='ebook')"
rm *.aux *.bbl *.blg *.log *.out *.toc *.dvi

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/
