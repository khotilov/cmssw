How to make PDF files from AI files for inclusion into documents:
 - Open file in Illustrator CS3
 - File / Save as...
    - Format: Adobe PDF (pdf), Save
    - Preset: Illustrator default (PDF 1.5 / Acrobat 6)
    - General:
       - Disable "Preserve Illustrator Editing Capabilities"
       - Enable "Embed Page Thumbnails"
       - Enable "Optimize for Fast Page View"
       - Disable "View PDF after Saving"
       - Disable "Create Acrobat Layers from Top-Level Layers"
    - Compression:
       - Everything set to "Do Not Downsample"
       - Enable "Compress Text and Line Art"
    - Marks and Bleeds: Disable everything, Marks Offset 6pt, Bleeds all at 0pt
    - Output: No color conversion, do not include profiles, PDF/X disabled
    - Advanced: Subset fonts when <100% used (= subset everything)
    - Security: No passwords or permissions set
    - Save PDF

How to convert Excel charts to PDF:
  - Open file in Office 2008 on a Mac
  - Select chart you want to save
  - Right click, Save as Picture...
  - Format: PDF, Select name, Save

  Note that gradients and sometimes also shadows produce
  PDF which can only be shown in Mac's Preview.  Always
  check the PDF displays correctly in Acrobat Reader.

How to convert PDF slide to a high-resolution transparent PNG:
  - Open PowerPoint 2008 on a Mac
  - PowerPoint / Preferences... / Save / Save slides as graphics files:
     - Save current slide only
     - 1200 dpi (or other desired resolution)
  - Optional: if this is an illustration-only slide, set page size in
    File / Page Setup such that the page just about wraps the graphics
  - Go to the slide you want to save
  - File / Save as Pictures..., Format: PNG, Select name, Save
  - Knock out the background from the PNG file in Photoshop:
     - Open file in Photoshop
     - Crop image if necessary
     - Select magic wand tool, tolerance zero, select contiguous
     - Click on background, and all other background areas such as
       closed loops in text
     - Select / Modify / Contract... by 2 pixels
     - Edit / Cut
     - Save

How to convert PowerPoint, Excel graphics to AI:
  - Open document in Office 2008 on a Mac
  - Select desired graphics elements
  - Edit / Copy
  - Open existing or new document in Illustrator CS3
  - Edit / Paste

  You may need to edit graphics. Text frequently gets broken down to typeset
  sequences, and needs to be merged back to full text boxes, including multi-
  line text.  Replace shadows and gradients with Illustrator's native ones,
  the ones from Office will import as special objects and will look bad.  You
  need to re-group and re-align everything. You get better quality PDFs if
  you first convert the Office graphics to Illustrator, save as PDF in AI.

How to extract original bitmap images from PDF, AI files without ripping:
  - Open PDF or AI document in Photoshop CS3
  - Select "Images" in the Import PDF dialog
  - Navigate to the image you want to extract
  - Edit / Convert to Profile... / sRGB (if necessary)
  - File / Save for Web & Devices..., Preset: PNG-24, Enable Transparency,
    Disable Interlaced, Save

Producing good quality PDF files from LaTeX using teTeX:
  - Create a temporary subdirectory, for example "Figures"
      THIS IS TEMPORARY, DO NOT COMMIT THESE FILES AS ORIGINALS!
      ORIGINALS SHOULD BE KEPT IN "Illustrations".

  - Create "safe" symlinked names and .bb files for the graphics:
    (cd Figures && for f in .../Illustrations/*.{pdf,png}; do
       dest=$(basename $f | tr ' ' _); [ -f "$dest" ] || ln -s $f $dest;
       case $f in *.png ) ebb $dest ;; esac
     done)

  - Verify all your PDF imports are at most PDF version 1.5 and that
    all image content uses safe colour profiles such as sRGB. If you
    followed the instructions above to produce files, you should be safe.

  - Use latex + dvipdfmx for optimal results. Include files from the "Figures"
    symlink using the underscore-based safe names, without file extension
    (file names with spaces are troublesome).  Use something like this in
    your .tex document:

     \usepackage[dvipdfm]{graphicx}
     \usepackage[dvipdfm]{hyperref}
     \usepackage{subfigure}
     \usepackage{mediabb}

     \begin{figure}[!tbp]
     \begin{center}
       \includegraphics[width=.9\textwidth]{Figures/DQM_end_to_end}
     \end{center}
     \caption{\label{fig:endtoend}Da picture.}
     \end{figure}

   - Process file with "latex", "bibtex", "makeindex", etc. as usual
   - Produce PDF with "dvipdfmx -V 5  -o $NAME.pdf $NAME.dvi"
   - Delete all unnecessary temporary files; keep only the .tex and .pdf
