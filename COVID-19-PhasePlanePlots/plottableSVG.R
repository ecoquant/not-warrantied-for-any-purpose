# plottableSVG.R
# Jan Galkowski, 27 June 2018
# Last changed 26th October 2020.

# via Inkscape

constructFilenameFrom<- function(root="pfx", suffix=".svg")
{
  stamp<- Sys.time()
  stamp<- gsub(x=stamp, pattern="(-|:)", replacement="", fixed=FALSE)
  stamp<- gsub(x=stamp, pattern=" ", replacement="-", fixed=TRUE)
  E<- round(100*proc.time()["elapsed"])
  return(list(full=sprintf("%s-%s-%s%s", root, stamp, E, suffix), stem=sprintf("%s-%s-%s", root, stamp, E)))
}

openSVG<- function(root="DefaultSVGStem", width=11, height=8, pointsize=12, antialias="gray", onefile=FALSE)
{
  fnx<- constructFilenameFrom(root=root, suffix=".svg")
  # (antialias="gray" is also an option below)
  svg(filename=fnx$full, width=width, height=height, pointsize=pointsize, onefile=onefile, family="mono", antialias=antialias)
  return(fnx)
}

closeSVG<- function(fnx)
{
  stopifnot(2 == length(fnx))
  stopifnot( all( c("full", "stem") %in% names(fnx) ) )
  dev.off()
  pathPfx<- getwd()
  shQCall<- shQuote("c:/Program Files/Inkscape/bin/inkscape")
  runResult<- system(command=sprintf("%s --export-type=\"pdf\" \"%s/%s\"", shQCall, pathPfx, fnx$full), 
#                     sprintf("%s.pdf", fnx$stem)), 
                 wait=TRUE, ignore.stdout=FALSE)
  unlink( sprintf("%s.pdf_tex", fnx$stem) )
  unlink( sprintf("%s.svg", fnx$stem) )
  return(runResult)
}

# Used as:
#\begin{figure}
#\centerline{\includegraphics[scale=0.75]{DefaultSVGStem-20180627-194106-41285.pdf}}
#\caption{\large TEST SVG}
#\label{fi:TestSVG}
#\end{figure}

