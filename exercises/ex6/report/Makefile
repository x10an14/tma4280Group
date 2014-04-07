
ex6.pdf: ex6.tex
	pdflatex ex6.tex
	bibtex ex6.aux
	pdflatex ex6.tex
	pdflatex ex6.tex #Two runs for references

clean:
	rm -f *.aux *.log *.pdf *.loc *.lot *.toc *.lof *-plot.eps *-plot.tex

all: clean ex6.pdf

clean:
	rm -f *.aux *.log *.pdf *.loc *.lot *.toc *.lof

