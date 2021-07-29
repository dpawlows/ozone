#!/bin/bash

#Delete folders created during last run of runMultipleEOM.sh
echo "Deleting folders from last set of runs"
rm -r runMultipleEOMdata/*

#Start log file to track data from all runs for plotting
echo '#The current order (by line) of the data is: (1) Pressure [P/PEarth], (2) Star Temp [K], (3) atmosphere temperature [K], (4) maximum ozone abundance [molecules/cm^3], (5) altitude of max abundance [km], (6) total column (molecules/cm^3)' >> runMultiplelog.txt
  #Start from the grep on line 34, then list all data you want from input*.inp

#Loop for input files you want to run in EOM
for inputFile in inputMultipleFiles/*.inp
do
  #Create a variable for each input file
  #echo $inputFile              #to print the inputFile variable
  temp="${inputFile##*/}"       #strip 'inputMultipleFiles/' from inputFile variable
  newdir="${temp%.*}"           #strip '.inp' from inputFile variable

  #Create a directory for iterations input file and directory to copy data after run is completed
  mkdir runMultipleEOMdata/$newdir
  mkdir runMultipleEOMdata/$newdir/data

  #Delete input files currently in input folder, Copy current input file for EOM to input folder (needed to run) and iteration folder
  rm input/*.inp
  cp $inputFile runMultipleEOMdata/$newdir/input.inp
  cp $inputFile input/input.inp

  #Run EOM, copy data to iteration data folder
  echo 'EOM is currently running with' $newdir'.inp as its input file'
  eom.py input > eomlog.txt

  #Gather data from input file and eom run
  grep '(Pressure Ratio)' $inputFile | sed 's/(.*)//' >> runMultiplelog.txt #Search for (P/PEarth) in $inputFile then remove any text between parenthesis in the line that matches (P/PEarth)
  grep '(TStar)' $inputFile | sed 's/(.*)//' >> runMultiplelog.txt
  #grep '(AU)' $inputFile | sed 's/(.*)//' >> runMultiplelog.txt
  #grep 'Value from input file you want to keep' $inputFile | sed 's/(.*)//' >> runMultiplelog.txt
    #Use the above line to pull which ever data you are interested in; adjust line 8 to reflect data you are keeping track of

  grep 'Atmosphere Temperature: ' eomlog.txt | tail -1 | sed 's/^.*: //' >> runMultiplelog.txt

  #Run (altered) plotTIme and plotProfile, Copy plots to iteration folder, append data to runMulitplelog.txt
  python plotprofileMultiple.py > profilelog.txt
  grep 'Max ozone: ' profilelog.txt | sed 's/^.*: //' >> runMultiplelog.txt #Search for Max ozone: in profilelog.txt substitue (s) from the beginning of the line (^) with any number of characters and any characters (.*) up to when the match ends (:) with the empty string; slashes are delimiters
  grep 'Max altitude: ' profilelog.txt | sed 's/^.*: //' >> runMultiplelog.txt
  grep 'Total column: ' profilelog.txt | sed 's/^.*: //' >> runMultiplelog.txt

  python plotTimeMultiple.py
  mv data/*.png runMultipleEOMdata/$newdir
  mv data/*.dat runMultipleEOMdata/$newdir/data

done

mv runMultiplelog.txt runMultipleEOMdata/runMultiplelog.txt
rm eomlog.txt profilelog.txt

python plotDataMultiple.py

mv *.png runMultipleEOMdata/
mv sortedData.txt runMultipleEOMdata/inputLog.txt
echo 'All data from this set of runs is available in the runMultipleEOMdata folder.'
