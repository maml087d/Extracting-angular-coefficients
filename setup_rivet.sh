if ! [[ -z $(screen -ls | grep -o rivetsetup1) ]]

then 
	screen -r rivetsetup1

else

	screen -d -m -L -Logfile /data/horse/ws/maml087d-workspace/logger.txt -S rivetsetup1
	screen -S rivetsetup1 -p 0 -X stuff "setupATLAS -c el9
	"

	screen -S rivetsetup1 -p 0 -X stuff "asetup 23.6.26,AthGeneration
	"

	screen -S rivetsetup1 -p 0 -X stuff "source setupRivet
	"

	screen -r rivetsetup1
fi
