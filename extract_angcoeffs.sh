#! /bin/bash
# extracting angular coeffs
# writing the coeff file.txt and extracting data histograms from yoda to json in outpath/histo_nocuts // outpath/histo_cuts

set -euo pipefail

# settings
cuts=true
fullrun=false
yodaname=Zjj_inclusive_data_cuts_lab.yoda.gz

ori_runfile=/data/horse/ws/maml087d-workspace/rivet_anas/run_ana.py
analysis=/data/horse/ws/maml087d-workspace/rivet_anas/angcoeff_inclusive_labframe.cc
# logfile=/data/horse/ws/maml087d-workspace/logger.txt

outpath=/data/horse/ws/maml087d-workspace/rivet_anas/data/labframe/inclusive
# IMPORTANT: outpath/histo_cuts or outpath/histo_nocuts must exist

ptbinfile=/data/horse/ws/maml087d-workspace/rivet_anas/data/labframe/inclusive/ptbins.txt

# creating copy of the original run file for this run and creating logger
declare -i f=0
while [ -e "rivet_anas/run_ana_$f.py" ]; do
    f+=1
done
runfile=/data/horse/ws/maml087d-workspace/rivet_anas/run_ana_$f.py
logfile=/data/horse/ws/maml087d-workspace/logger$f.txt
echo "runfile: $runfile; log: $logfile"

cp $ori_runfile $runfile

touch $logfile
 > $logfile

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

# making screen session
name=rivet-ang
declare -i j=0
while ! [[ -z $(screen -ls | grep -o $name$j) ]]; do
	j+=1
done
screenname=$name$j

screen -d -m -L -Logfile $logfile -S $screenname
echo "made screen: $screenname"

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
screen -S $screenname -p 0 -X stuff "rm $runfile
"

# reading the output angcoeffs from yoda file
if $cuts;
then
	echo "with cuts"
	screen -S $screenname -p 0 -X stuff "python read_angcoeffs.py $yodafile $outpath/histo_cuts /ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/ $ptbinfile
	"
	screen -S $screenname -p 0 -X stuff "python read_binned_temps.py $yodafile $outpath/histo_cuts /ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/ $ptbinfile
	"
else
	echo "without cuts"
	screen -S $screenname -p 0 -X stuff "python read_angcoeffs.py $yodafile $outpath/histo_nocuts /ATLAS_2020_I1803608:CUTS=NO:TYPE=EW_ONLY/ $ptbinfile
	"
	screen -S $screenname -p 0 -X stuff "python read_binned_temps.py $yodafile $outpath/histo_nocuts /ATLAS_2020_I1803608:CUTS=NO:TYPE=EW_ONLY/ $ptbinfile
	"
# 	screen -S $screenname -p 0 -X stuff "python read_binned_temps.py $yodafile $outpath /ATLAS_2020_I1803608:CUTS=NO:TYPE=EW_ONLY/ $ptbinfile
# 	"
fi

screen -S $screenname -p 0 -X stuff "exit
"

# closing apptainer and job
if $fullrun
then
	screen -S $screenname -p 0 -X stuff "exit
	"
fi
