#!/bin/bash


### Script for convergence.


### Valid Grid size:
###
### {3 4 5 6 7 8 9 10 
###  11 12 13 14 15 17 19 21
###  23 25 27 29 33 37 41 45 
###  49 53 57 65 73 81 89 97 
###  105 113 129 145 161 177 193 209 
###  225 257 289 321 353 385 417 449 
###  513 577 641 705 769 833 897 1153 
###  1281 1409 1665 1793 2305 2817 3329};

            
order=6



echo "Runing Olliptic for convergence..."

NRes="41 45 49 53 57 65 73 81 161"


problem="Punctures"
testsystem=0           ### 1 for true 0 for false


directory="SBH"
box_size=100
levels=6

bhole="-NumBH 1 -bhmp1 1.0 -bhsz1 0.2"

boundary="-B Robin"
#boundary="-B Asymptotic -BA 0.0 -Bn 1"

#symmetry="-ds octant"


tols=1e-10
par="-nu2 20"
Dprint=0.01

exebin="./Olliptic.x"


##########

conv_dir="Conv"$directory"L"$box_size"O"$order"l"$levels
echo "set term wxt" > $conv_dir".gnuplot"

echo "output to: "$conv_dir


num=1
for N in `echo $NRes`
do


    
output=$directory"L"$box_size"R"$num"O"$order"l"$levels

stdout=$output".out"
stderr=$output".err"


echo "=============="
echo "starting $N..."
time $exebin -eq $problem $boundary $bhole -L $box_size -N $N -O $order -l $levels  -s $tols -Dp $Dprint -o $output $par $symmetry > $stdout 2> $stderr
mv $stdout $output
mv $stderr $output
echo "$N done"
echo "=============="

echo "N"$num"="$N >> $conv_dir".gnuplot"
h=$(echo "scale=9; $box_size/($N-1)" | bc)
echo "h"$num"=0"$h >> $conv_dir".gnuplot"




((num++))

done


((NumMax=num))


mkdir $conv_dir
mv $directory"L"$box_size*"O"$order"l"$levels $conv_dir
mv $conv_dir".gnuplot" $conv_dir
cd $conv_dir


num=1


((numM=NumMax-1))
while [ "$num" -lt "$numM" ]
do


((nump1=num+1))


if [ "$num" == "1" ]; then
output1=$directory"L"$box_size"R"$num"O"$order"l"$levels"/u.X"
else
output1="u"$num".X"
fi
output2=$directory"L"$box_size"R"$nump1"O"$order"l"$levels"/u.X"



tmp="u"$nump1".X"

join $output1 $output2 > $tmp

((num++))
done

joinfile="u"$conv_dir".X"

mv $tmp $joinfile

echo "p=$order ">> $conv_dir".gnuplot"
echo "set title sprintf(\"|*u-u^h|\\n___________\\n|*h^%1.3f-h^%1.3f|\",p,p)" >> $conv_dir".gnuplot"
echo "plot \\">> $conv_dir".gnuplot"

if [ "$testsystem" == "1" ]; then
N=2
((NNmax=2*(NumMax-1)))
dN=2
else
N=2
((NNmax=NumMax))
dN=1
fi
((numbase=(NumMax-1) ))

num=1
while [ $N -lt $NNmax ]
do

echo -n "\"$joinfile\" u 1:(abs(\$$NNmax-\$$N)/abs(h$numbase**p-h$num**p)) w l t sprintf(\"h=%f,N=%1.0f\",h$num,N$num)">> $conv_dir".gnuplot"
((num++))
((N+=dN))

if [ "$N" -lt $NNmax ]; then
echo ",\\">> $conv_dir".gnuplot"
else
echo " " >> $conv_dir".gnuplot"
fi


done

echo "pause -1" >> $conv_dir".gnuplot"
echo "set term postscript eps color" >> $conv_dir".gnuplot"
echo "set output \"SelfConv.eps\" " >> $conv_dir".gnuplot"
echo "replot " >> $conv_dir".gnuplot"



echo "f(x)=a+b*x" >> "ConvS.gnuplot"
echo "set fit logfile \"CSfit.log\" " >> "ConvS.gnuplot"
echo "set fit errorvariables" >> "ConvS.gnuplot"
echo "fit f(x) \"uL1.dat\" via a,b" >> "ConvS.gnuplot"
echo "set ylabel \"Log(*E(h)) (adim)\" " >> "ConvS.gnuplot"
echo "set xlabel \"Log(h) (adim)\" " >> "ConvS.gnuplot"
echo "set key bottom" >> "ConvS.gnuplot"
echo "set title sprintf(\"*E(h)~C h^ p :\n Log(C)= %f +/- %f \n p = %f +/- %f\",a,a_err,b,b_err)" >> "ConvS.gnuplot"
echo "plot \"uL1.dat\" t \"Olliptic\", a+b*x t sprintf(\"%f + %f x\",a,b) " >> "ConvS.gnuplot"
echo "pause -1" >> "ConvS.gnuplot"
echo "set term postscript eps color" >> "ConvS.gnuplot"
echo "set output \"ConvS.eps\" " >> "ConvS.gnuplot"
echo "replot " >> "ConvS.gnuplot"




if [ "$testsystem" == "1" ]; then

echo "set term wxt">> $conv_dir".gnuplot"
echo "set title sprintf(\"|u-u^h|\\n___________\\nh^%1.3f\",p)" >> $conv_dir".gnuplot"
echo "plot \\">> $conv_dir".gnuplot"


col=3
num=1
for N in `echo $NRes`

do

output1=$directory"L"$box_size"R"$num"O"$order"l"$levels
echo -n "\"$output1/u.X\" u 1:(abs(\$$col)/h$num**p) w l t sprintf(\"h=%f,N=%1.0f\",h$num,N$num)">> $conv_dir".gnuplot"
((num++))

if [ "$num" -lt $NumMax ]; then
echo ",\\">> $conv_dir".gnuplot"
else
echo " " >> $conv_dir".gnuplot"
fi

done

echo "pause -1">> $conv_dir".gnuplot"
echo "set term postscript eps color" >> $conv_dir".gnuplot"
echo "set output \"ErrorConv.eps\" " >> $conv_dir".gnuplot"
echo "replot " >> $conv_dir".gnuplot"


echo "f(x)=a+b*x" >> "ConvA.gnuplot"
echo "set fit logfile \"CAfit.log\" " >> "ConvA.gnuplot"
echo "set fit errorvariables" >> "ConvA.gnuplot"
echo "fit f(x) \"vL1.dat\" via a,b" >> "ConvA.gnuplot"
echo "set ylabel \"Log(E(h)) (adim)\" " >> "ConvA.gnuplot"
echo "set xlabel \"Log(h) (adim)\" " >> "ConvA.gnuplot"
echo "set key bottom" >> "ConvA.gnuplot"
echo "set title sprintf(\"E(h)~C h^ p :\n Log(C)= %f +/- %f \n p = %f +/- %f\",a,a_err,b,b_err)" >> "ConvA.gnuplot"
echo "plot \"vL1.dat\" t \"Olliptic\", a+b*x t sprintf(\"%f + %f x\",a,b) " >> "ConvA.gnuplot"
echo "pause -1" >> "ConvA.gnuplot"
echo "set term postscript eps color" >> "ConvA.gnuplot"
echo "set output \"ConvA.eps\" " >> "ConvA.gnuplot"
echo "replot " >> "ConvA.gnuplot"




fi




echo "...all done"


