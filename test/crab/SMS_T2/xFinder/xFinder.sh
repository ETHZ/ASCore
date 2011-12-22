#!/bin/bash

##you may want to change the next line.
NEvents=123


collectionnum=0
#masses
mglu[$NEvents]=0
mlsp[$NEvents]=0
mchi[$NEvents]=0

#number of entries
nmglu[$NEvents]=0
nmlsp[$NEvents]=0
nmchi[$NEvents]=0

chatty=0

function analyzeline {
    if [[ "$1" == *ParticleListDrawer* ]]; then
      let 'collectionnum = collectionnum +1'
      if [[ $chatty -gt 0 ]]; then 
	echo "Found line with ParticleListDrawer; ($1) incremented collectionnum to $collectionnum"
      fi
    fi
    if [[ $1 == *1000021* ]]; then 
        result=`echo $1 |  awk '{ print $24 }'`
        Nresult=1
        oldresult=${mglu[collectionnum]}
        oldNresult=${nmglu[collectionnum]}
        if [ -n "${mglu[collectionnum]}" ]; then
        	newresult=`echo $result + $oldresult | bc`
        	newNresult=`echo $Nresult + $oldNresult | bc`
	else 
        	newresult=$result
        	newNresult=$Nresult
	fi
        mglu[$collectionnum]=$newresult
        nmglu[$collectionnum]=$newNresult
        if [[ $chatty -gt 0 ]]; then 
		echo "Found gluino mass. Incremented to ${mglu[collectionnum]} with  ${nmglu[collectionnum]}"
        fi
    fi
    if [[ $1 == *1000023* ]]; then 
        result=`echo $1 |  awk '{ print $24 }'`
        Nresult=1
        oldresult=${mchi[collectionnum]}
        oldNresult=${nmchi[collectionnum]}
        if [ -n "${mchi[collectionnum]}" ]; then
        	newresult=`echo $result + $oldresult | bc`
        	newNresult=`echo $Nresult + $oldNresult | bc`
	else 
        	newresult=$result
        	newNresult=$Nresult
	fi
        mchi[$collectionnum]=$newresult
        nmchi[$collectionnum]=$newNresult
        if [[ $chatty -gt 0 ]]; then 
		echo "Found chi mass. Incremented to ${mchi[collectionnum]} with  ${nmchi[collectionnum]}"
        fi
    fi
    if [[ $1 == *1000022* ]]; then 
        result=`echo $1 |  awk '{ print $24 }'`
        Nresult=1
        oldresult=${mlsp[collectionnum]}
        oldNresult=${nmlsp[collectionnum]}
        if [ -n "${mlsp[collectionnum]}" ]; then
        	newresult=`echo $result + $oldresult | bc`
        	newNresult=`echo $Nresult + $oldNresult | bc`
	else 
        	newresult=$result
        	newNresult=$Nresult
	fi
        mlsp[$collectionnum]=$newresult
        nmlsp[$collectionnum]=$newNresult
        if [[ $chatty -gt 0 ]]; then 
		echo "Found lsp mass. Incremented to ${mlsp[collectionnum]} with  ${nmlsp[collectionnum]}"
        fi
    fi
}



cipher=`date +"%Y_%m_%d_%H_%M_%S"`
echo "Fetching information about one thousand events...."
outfile="Analyzing_Particles_${cipher}.txt"
cmsRun ParticleListDrawer.py nEvents=$NEvents > $outfile
echo "Information fetched. Now let's process it."


while read line
do
analyzeline "$line"
done < "$outfile" 

echo "Analyzed all events and filled mglu,mchi, and mlsp arrays. Will compute x values now."
runover="1..$collectionnum"

sumx=0
nx=0

for (( i=1; i <= $collectionnum; i++ ))
do
	massglu=${mglu[i]}
	masschi=${mchi[i]}
	masslsp=${mlsp[i]}
	nglu=${nmglu[i]}
	nchi=${nmchi[i]}
	nlsp=${nmlsp[i]}
	massglu=$(echo "$massglu / $nglu" |bc -l)
	masschi=$(echo "$masschi / $nchi" |bc -l)
	masslsp=$(echo "$masslsp / $nlsp" |bc -l)
	x=$(echo "($masschi-$masslsp) / ($massglu-$masslsp)" |bc -l)
	echo "For step $i we have    G:$massglu  [$nglu], C:$masschi   [$nchi], L:$masslsp   [$nlsp] ---> x=$x"
	sumx=`echo $sumx + $x | bc`
	nx=`echo $nx + 1 | bc`
done

finalx=$(echo "$sumx / $nx" | bc -l)
echo "FINAL RESULT : x = $finalx , based on the first $nx events."
rm $outfile
exit