#
# makefile
# Field G. Van Zee
#

DOC := libflame

PDFLATEX  := /lusr/bin/tex/pdflatex
BIBTEX    := /lusr/bin/tex/bibtex
MAKEINDEX := /lusr/bin/tex/makeindex

all: pdf

latex:
	latex $(DOC).tex
	bibtex $(DOC)
	latex $(DOC).tex
	latex $(DOC).tex

ps: latex
	dvips -f $(DOC).dvi -t letter > $(DOC).ps

pdf:
	$(PDFLATEX)  $(DOC).tex
	$(BIBTEX)    $(DOC)
	$(MAKEINDEX) $(DOC)
	$(PDFLATEX)  $(DOC).tex
	$(PDFLATEX)  $(DOC).tex

clean:
	- rm -f $(DOC).aux $(DOC).idx $(DOC).log $(DOC).blg $(DOC).bbl $(DOC).toc
	- rm -f $(DOC).dvi $(DOC).pdf $(DOC).ps $(DOC).ind
	- rm -f make.err
	- find . -name "*~" | xargs rm -f
