all: vsn likelihoodcomputations convergence2 clean

vsn: vsn.tex
	texi2dvi --pdf vsn.tex

likelihoodcomputations: likelihoodcomputations.tex
	texi2dvi --pdf likelihoodcomputations.tex

convergence2: convergence2.tex
	cp -p ../scripts/convergence2.pdf .

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
	-$(RM) -f vsn-* likelihoodcomputations-* Rplots.pdf



