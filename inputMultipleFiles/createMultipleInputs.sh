#echo "How many copies would you like to make?"
#read copies

let i=1
#let end=$copies+1

temperature=(2500 3500 4500 5500)
pressure=(0.1 1 5 10 25 50)

#while ((i < $end))
#do
  #iFile=$(printf "input%.03d.inp" $i)
  #echo $iFile
  #cp /Users/blaycock/Research/ozone/inputFiles/input000.inp /Users/blaycock/Research/ozone/inputMultipleFiles/$iFile
  #sed -i '' "s/.*(time step,/$tstep (time step, /g" "$iFile"
  #let i++
  #tstep=$(($tstep + 150))
#done

for temps in ${!temperature[@]}
do
  for press in ${!pressure[@]}
  do
    iFile=$(printf "input%.03d.inp" $i)
    cp /Users/blaycock/Research/ozone/inputFiles/input000.inp /Users/blaycock/Research/ozone/inputMultipleFiles/$iFile
    sed -i '' "s/.*(Tstar)/${temperature[$temps]} (TStar)/g" "$iFile"
    sed -i '' "s/.*(Pressure Ratio)/${pressure[$press]} (Pressure Ratio)/g" "$iFile"
    let i++
    echo $iFile
  done
done
