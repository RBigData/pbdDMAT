#!/bin/sh

PKGVER=`grep "Version:" ../DESCRIPTION | sed -e "s/Version: //"`
sed -i -e "s/myversion{.*}/myversion{${PKGVER}}/" pbdDMAT-guide.Rnw


rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pbdDMAT-guide.Rnw
bibtex pbdDMAT-guide
#makeindex pbdDMAT-guide
#pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
Rscript -e "tools::compactPDF('pbdDMAT-guide.pdf', gs_quality='ebook')"
rm *.aux *.bbl *.blg *.log *.out *.toc *.dvi

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/
