.PHONY: all build deploy

all: build

build: ../dist/gh-pages/spEnergiesQD/figures.pdf
	cd Manuscript && $(MAKE)

deploy: ../dist/gh-pages/.git/config build
	cd ../dist/gh-pages && \
	git add -A && \
	git commit --amend -q -m Autogenerated && \
	git push -f origin master:gh-pages

../dist/gh-pages/.git/config:
	mkdir -p ../dist/gh-pages
	url=`git remote -v | grep origin | awk '{ printf "%s", $$2; exit }'` && \
	cd ../dist/gh-pages && \
	git init && \
	git config user.name Bot && \
	git config user.email "<>" && \
	git commit -m _ --allow-empty && \
	git remote add origin "$$url"

../dist/gh-pages/spEnergiesQD/figures.pdf: figures.pdf
	mkdir -p $(@D)
	cp figures.pdf $@

.SUFFIXES: .md .pdf .svg .tex
.svg.pdf:
	svg2pdf $<

.md.tex:
	pandoc -s -V fontsize=12pt -o $@ $<

.tex.pdf:
	latexmk -interaction=nonstopmode -pdf $*

figures.dep: figures.md gen-deps
	./gen-deps figures

-include figures.dep

.SECONDARY:
