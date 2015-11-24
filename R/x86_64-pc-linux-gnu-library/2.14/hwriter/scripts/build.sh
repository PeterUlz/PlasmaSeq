#!/bin/sh

rm -rf hwriter.Rcheck .Rhistory ..Rcheck example-hwriter.html inst/doc/hwriter.pdf hwriter*.tar.gz article

R CMD build .
R CMD check hwriter*.tar.gz
rm -rf hwriter.Rcheck

R CMD INSTALL . 
echo "library(hwriter) ; example(hwriter) ; file.copy(file.path(tempdir(),'example-hwriter.html'),'.',overwrite=TRUE)" | R --no-save --vanilla

scp -P 6422 hwriter_1.2.tar.gz gpau@localhost:~/public_html/hwriter
scp -P 6422 example-hwriter.html gpau@localhost:~/public_html/hwriter/index.html

rm -rf hwriter.Rcheck .Rhistory ..Rcheck inst/doc/hwriter.pdf
