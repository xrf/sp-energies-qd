all: ../../dist/gh-pages/spEnergiesQD/Manuscript/paper.pdf

clean:
	rm -fr *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out *Notes.bib paper.pdf

.PHONY: all clean

../../dist/gh-pages/spEnergiesQD/Manuscript/paper.pdf: paper.pdf
	mkdir -p $(@D)
	cp paper.pdf $@

paper.pdf: paper.bib paper.tex head_plain.tex head_revtex.tex
	../mklatex paper
