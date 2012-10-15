#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex pbdDMAT-guide.Rnw
bibtex pbdDMAT-guide
pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
pdflatex pbdDMAT-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc *.dvi
