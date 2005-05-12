
*** Skeleton of the note on Raw Data Formats ***

This directory (DataFormats/domain/doc/RawDataFormats) contains the
master files for the note describing the primary raw data formats.

The master .tex file is DataFormats.tex. This can be compiled with
both latex and pdflatex (recommended).


The section describing the format for each detector is located in a
separate file, named RawDataFormat.tex, within the doc/tex directory
of the relative detector package:

DataFormats/SiPixelRawData/doc/tex/RawDataFormat.tex
DataFormats/SiStripRawData/doc/tex/RawDataFormat.tex
DataFormats/EcalRawData/doc/tex/RawDataFormat.tex
DataFormats/HcalRawData/doc/tex/RawDataFormat.tex
DataFormats/DTRawData/doc/tex/RawDataFormat.tex
DataFormats/CSCRawData/doc/tex/RawDataFormat.tex
DataFormats/RPCRawData/doc/tex/RawDataFormat.tex

All of these must be present to be able to compile DataFormats.tex.

To help preparing tables with data format specifications, a set of
macros is available. An example of their usage is included at the end
of the note (source file Examples.tex). This is of course just a base,
it is not meant to cover all possible needs.

In case tables with data formats are already existing in a different
format (e.g. .doc) and you don't want to retype them, they can be
included as eps/pdf objects. Please contact us if you need help in
doing this.


For questions, feedback, requests, help: 
nicola.amapane@cern.ch
emilio.meschi@cern.ch
