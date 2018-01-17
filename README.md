# Addition and removal energies of circular quantum dots

  - [Latest draft](https://xrf.github.io/sp-energies-qd)
  - [arXiv:1707.00229](https://arxiv.org/abs/1707.00229)
  - [doi:10.1063/1.4995615](https://doi.org/10.1063/1.4995615)

## Building

Dependencies:

- [GNU Make](https://www.gnu.org/software/make/)
- [LaTeX](https://www.latex-project.org) (must have various packages as well as `latexmk`)
- [Python](https://www.python.org) (for Makefile dependency generation and other scripts)

Optional dependencies (needed for building `figures.pdf`):

- [Inkscape](https://inkscape.org) (for SVG to PDF conversion)
- [Pandoc](https://pandoc.org) (for converting Markdown to LaTeX)

To build just the paper, simply run:

~~~sh
make
~~~

The PDF will be output to `Manuscript/paper.pdf`.

To build the notes (`figures.pdf`), run:

~~~sh
make notes
~~~
