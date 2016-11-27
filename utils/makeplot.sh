#!/bin/bash
E_BADARGS=65

#sizeX="5.2in"
#sizeY="5.2in"

sizeX="5.2in"
sizeY="3.7in"

if [ -z "$1" ]
then
  echo "Usage: `basename $0` filename.gplot"
  exit $E_BADARGS
fi  

if [ ! -f "$1" ] 
then
 echo "File $1 not found."
  exit $E_BADARGS
fi  

extname=`ls $1 | cut -f 2 -d "."`
echo "File $1 found."

if [ ! "$extname" == "gplot" ]  
then
  echo "Usage: `basename $0` filename.gplot"
  exit $E_BADARGS          # Exit and explain usage.
fi  

echo "making plot: gnuplot $1"

if ! `gnuplot $1`; then echo "command failed"; exit 1; fi 

echo "...done"


filename=`ls $1 | cut -f 1 -d "."`

texfile=$filename".tex"
epsfile=$filename".eps"


if [ ! -f $texfile ] 
then
 echo "File $texfile not found."
  exit $E_BADARGS
fi  

if [ ! -f $epsfile ] 
then
 echo "File $epsfile not found."
  exit $E_BADARGS
fi  

echo "Files $epsfile and $texfile found."

echo "Making pdf..."


tmpfile=$filename"_tmp.tex"
dvifile=$filename"_tmp.dvi"
pdffile=$filename"_tmp.pdf"

echo $filename


echo "\documentclass[10pt]{article}">> $tmpfile


echo "\usepackage[margin=0in, paperwidth=$sizeX, paperheight=$sizeY]{geometry}">> $tmpfile
echo "\usepackage{nopageno}">> $tmpfile
echo "\usepackage{txfonts}">> $tmpfile
echo "\usepackage[usenames]{color}">> $tmpfile

echo "\usepackage{graphicx}">> $tmpfile
echo "\usepackage{latexsym}">> $tmpfile
echo "\usepackage[draft=false]{hyperref}">> $tmpfile

echo "\definecolor{rgb1}{RGB}{65,105,225}">> $tmpfile
echo "\definecolor{rgb2}{RGB}{178, 34, 34}">> $tmpfile
echo "\definecolor{rgb3}{RGB}{46,139, 87}">> $tmpfile

echo "\begin{document}">> $tmpfile

echo "\input{ $texfile }">> $tmpfile 

echo "\end{document}">> $tmpfile


texi2dvi $tmpfile > $filename".log"
dvipdf $dvifile

echo "Cleaning..."

mv $pdffile $filename".pdf"
rm *tmp*
rm $filename".log"
rm $filename".eps"
rm $filename".tex"

echo "...all done."

