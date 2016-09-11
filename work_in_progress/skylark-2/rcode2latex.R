library(stringr)
rm(list=ls())

latex      = "^\\s*# (.*)"
section       = "^\\s*#\\s?#+\\s+(.*)\\s+#{3,}$"
subsection    = "^\\s*#\\s?=+\\s+(.*)\\s+={3,}$"
subsubsection = "^\\s*#\\s?-+\\s+(.*)\\s+-{3,}$"
# section       = "^\\s*#1\\s+(.*)"
# subsection    = "^\\s*#2\\s+(.*)"
# subsubsection = "^\\s*#3\\s+(.*)"
ignore     = "^\\s*(#[-=]+)|(##.*)|(#'.*)"
ignore2    = "^\\s*##.*"
empty      = "^\\s*$"
mixed      = "^(.*[[:print:]].*)#(.+)"
quit       = "^\\s*#stop.*"

chunks = c("workhorse", "gof", "post","wald","index","overall")
lines = character(0)
for (s in chunks) {
  name = sprintf("trim_%s.R", s)
  new.lines <- trimws(readLines(name), "right")
  lines <- c(lines, new.lines)
}

texname <- "trim_code.tex"
tex = file(texname, "wt")

header="\\documentclass[a4paper]{article}
\\usepackage{amsmath}
\\usepackage{mathtools} % defines pmatrix*
\\usepackage{verbatim}
\\usepackage{parskip}
\\usepackage[margin=72pt]{geometry}
\\usepackage{color}
\\begin{document}
\\let\\Sum=\\sum
\\newcommand{\\diag}[1]{\\operatorname{diag}(#1)}
\\newcommand{\\var}[1]{\\operatorname{var}(#1)}
\\newcommand{\\cov}[1]{\\operatorname{cov}(#1)}
\\newcommand{\\se}[1]{\\operatorname{S.E.}(#1)}
\\newcommand{\\Mu}{{\\mu_{+}}}
\\title{TRIMR --- TRIM in R}
\\author{Patrick Bogaart, Mark van der Loo, Jeroen Pannekoek, Arco van Strien}
\\date{\\today}
\\maketitle
\\tableofcontents"

cat(header, "\n", file=tex)

old.mode    = '-'
cur.mode    = 'T'
special=""

set.mode <- function(new.mode, space=TRUE) {
  modes = paste0(cur.mode, new.mode)
  if (modes=="RR") cat("\\newline\n", file=tex)
  if (modes=="RT") cat("\\par\n", file=tex)
  if (modes=="RE") cat("\\par\n", file=tex)
  if (modes=="TT") cat("\n", file=tex)
  if (modes=="TR") cat("\\par\n", file=tex)
  if (modes=="TE") cat("\\par\n", file=tex)
  if (modes=="EE") cat("\n", file=tex)
  if (modes=="ET") cat("\n", file=tex)
  if (modes=="ER") cat("\n", file=tex)
  cur.mode <<- new.mode
}

go.R = function() {
  if (mode!='R') {
    cat("\\begin{verbatim}\n", file=tex)
    mode <<- "R"
  }
}

go.tex = function() {
  if (mode!='T') {
    cat("\\end{verbatim}\n", file=tex)
    mode <<- "T"
  }
}

for (line in lines) {

  if (grepl(quit, line)) {
    break
  }
  else if (grepl(section, line)) {
    set.mode('T', FALSE)
    txt = str_match(line, section)[1,2]
    cat(sprintf("\n\n\\section{%s}", txt), file=tex)
  }
  else if (grepl(subsection, line)) {
    set.mode('T', FALSE)
    txt = str_match(line, subsection)[1,2]
    cat(sprintf("\n\n\\subsection{%s}", txt), file=tex)
  }
  else if (grepl(subsubsection, line)) {
    set.mode('T', FALSE)
    txt = str_match(line, subsubsection)[1,2]
    cat(sprintf("\n\n\\subsubsection{%s}", txt), file=tex)
  }
  else if (grepl(latex, line)) { # R comment -> latex
    set.mode('T')
    txt = str_match(line, latex)[1,2]
    cat(txt, file=tex)
  }
  else if (grepl(ignore,line)) {
    # pass
  }
  else if (grepl(empty, line)) {
    set.mode('E')
    #cat("\\par\n", file=tex)
  }
  else if (grepl(mixed, line)) {
    set.mode('R')
    m = str_match(line, mixed)
    Rcode = trimws(m[1,2], "right")
    Lcode = trimws(m[1,3], "left")
    cat(sprintf("\\verb~> %s  ~{\\sffamily\\# %s}", Rcode, Lcode), file=tex)
  }
  else { # R code
    set.mode('R')
    cat(sprintf("\\verb~> %s~", line), file=tex)
    indent <- 1
    while (line[indent]==" ") indent <- indent+1
  }
}

trailer="\\end{document}\n"
cat(trailer, file=tex)
close(tex)

system(sprintf("pdflatex --interaction=nonstopmode %s", texname))
