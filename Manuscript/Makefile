../../dist/gh-pages/spEnergiesQD/Manuscript/paper.pdf: paper.pdf
	mkdir -p $(@D)
	cp paper.pdf $@

# TODO: stop ignoring errors
paper.pdf: paper.bib paper.tex
	latexmk -interaction=nonstopmode -pdf $* || :
