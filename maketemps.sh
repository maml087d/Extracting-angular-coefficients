#! /bin/bash
# makes templates

# settings
cuts=false
fullrun=true
yodaname=Zjj_ptbin_temps_nocuts.yoda.gz
runfile=/data/horse/ws/maml087d-workspace/rivet_anas/run_ana.py
analysis=/data/horse/ws/maml087d-workspace/rivet_anas/angcoeff.cc
logfile=/data/horse/ws/maml087d-workspace/logger.txt

sed -i "s/rivet.HistoFile = '.*'/rivet.HistoFile = '$yodaname'/" $runfile
# getting yoda file name from run_ana
yodafile=$(grep "^rivet.HistoFile" $runfile | sed "s/.*= '\(.*\)'/\1/")
echo "yodafile: $yodafile"

# setting the event number in run_ana.py
if $fullrun
then
	sed -i '1s/.*/theApp.EvtMax = -1/' $runfile
	echo "full run"
else
	sed -i '1s/.*/theApp.EvtMax = 10000/' $runfile
	echo "10000 events"
fi

# adjusting cut settings in run_ana
line=$(grep -m 1 -n "rivet.Analyses" $runfile | cut -d: -f1)
echo "line $line"

if $cuts
then
	sed -i "${line}s/.*/rivet.Analyses += [ '\ATLAS_2020_I1803608:TYPE=EW_ONLY:CUTS=YES\' ]/" $runfile
else
	sed -i "${line}s/.*/rivet.Analyses += [ '\ATLAS_2020_I1803608:TYPE=EW_ONLY:CUTS=NO\' ]/"  $runfile
fi

# setting up output directory for histograms
folder="tmp"
declare -i i=0
while [[ -d "/data/horse/ws/maml087d-workspace/templates/$folder$i" ]]; do
	i+=1
done

mkdir /data/horse/ws/maml087d-workspace/templates/$folder$i
outpath=/data/horse/ws/maml087d-workspace/templates/$folder$i
echo $outpath

# making screen session
name=rivet
declare -i j=0
while ! [[ -z $(screen -ls | grep -o $name$j) ]]; do
	j+=1
done 
screenname=$name$j
echo $screenname
screen -d -m -L -Logfile $logfile -S $screenname

echo "made screen"

if $fullrun
then
	screen -S $screenname -p 0 -X stuff "inter_long
	"
fi
# setting up atlas env in screen and running the $analysis analysis
screen -S $screenname -p 0 -X stuff "setupATLAS -c el9
"

screen -S $screenname -p 0 -X stuff "asetup 23.6.26,AthGeneration
"

screen -S $screenname -p 0 -X stuff "source setupRivet
"
screen -S $screenname -p 0 -X stuff "cd rivet_anas
" 
screen -S $screenname -p 0 -X stuff "rivet-build $analysis && athena $runfile
"

# reading the output histogram from yoda file
if $cuts;
then
	echo "with cuts"
	screen -S $screenname -p 0 -X stuff "python read_yodas_auto.py $yodafile  $outpath /ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/ data_2d
	"
	screen -S $screenname -p 0 -X stuff "python read_yodas_auto.py $yodafile  $outpath /ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/ 
	"
else
	echo "without cuts"
	screen -S $screenname -p 0 -X stuff "python read_yodas_auto.py $yodafile $outpath /ATLAS_2020_I1803608:CUTS=NO:TYPE=EW_ONLY/ data_2d
	"
	screen -S $screenname -p 0 -X stuff "python read_yodas_auto.py $yodafile $outpath /ATLAS_2020_I1803608:CUTS=NO:TYPE=EW_ONLY/
	"
fi

screen -S $screenname -p 0 -X stuff "exit
"

# closing apptainer and job
if $fullrun
then
	screen -S $screenname -p 0 -X stuff "exit
	" 
fi

