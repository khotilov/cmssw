#!/bin/sh


################# change here ###############
OUTFILE=CMS_IN_2009_666
TEXFILES="title.tex introduction.tex conclusion.tex acknow.tex bib.tex"
#############################################


echo "Merging files: $TEXFILES"
echo "               into the output file ${OUTFILE}.tex"

cat > ${OUTFILE}.tex <<EOF
\documentclass[colclass=cmspaper]{combine}
\usepackage{lineno}
\usepackage{amsfonts,amsmath,amssymb}
\usepackage[dvips]{graphicx}
\usepackage{bm}
\begin{document}
\begin{linenumbers}
\pagestyle{combine}
EOF

cat $TEXFILES >> ${OUTFILE}.tex

cat >> ${OUTFILE}.tex <<EOF
\end{linenumbers}
\end{document}
EOF

latex ${OUTFILE}.tex
bibtex ${OUTFILE}.aux
latex ${OUTFILE}.tex
latex ${OUTFILE}.tex
dvipdf ${OUTFILE}.dvi ${OUTFILE}.pdf
