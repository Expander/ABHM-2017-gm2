# http://www.gnu.org/software/make/manual/make.html

OUTFILENAME := talk.pdf
PLOTS       := 
TEXDIRS     := $(PLOTSDIR)
BIBTEX      := bibtex

.PHONY: all clean

all: $(OUTFILENAME)

$(PLOTS): plots/FlexibleEFTHiggs/*.dat plots/FlexibleEFTHiggs/MSSMtower-1L-nologs/*.dat plots/FlexibleEFTHiggs/*.sh plots/FlexibleEFTHiggs/*.gnuplot
	cd plots/FlexibleEFTHiggs && ./plot.sh

%.pdf: %.tex $(PLOTS)
	pdflatex $<
	cd Feynman && ./makeall.sh
	pdflatex $<

clean:
	rm -f *~ *.bak *.aux *.log *.toc *.bbl *.blg *.nav *.out *.snm *.backup
	rm -f plots/*.aux plots/*.log
	rm -f $(PLOTS)

distclean: clean
	rm -f $(OUTFILENAME)
