#!/bin/bash

if [ $# != 3 ];
then
  echo "USAGE: csta freq eruption_cycle "
  exit 
fi

ND=`pwd`
csta=$1
freq=$2
erup=$3
in_path=ERUPTION/${csta}/${freq}hz/${erup}
infile=${in_path}.txt
infold=${in_path}

if [ ! -d "Fig_ppm_single" ];then
	mkdir Fig_ppm_single/
fi
if [ ! -d "Fig_ppm_single/${csta}" ];then
	mkdir Fig_ppm_single/${csta}
fi
if [ ! -d "Fig_ppm_single/${csta}/${erup}" ];then
        mkdir Fig_ppm_single/${csta}/${erup}
fi
ofold=Fig_ppm_single/${csta}/${erup}

line=`cat $infile | wc -l`

eruption=`echo $erup | awk -F/ '{print$4}'`
#for ((i=1;i<=$line;i++))
for ((i=1;i<=$line;i+=5))
do
	rm ${csta}-???_Z?*.tmp
	stage=`awk 'NR==i{print $2}' i=${i} $infile`
	tt=`awk 'NR==i{print $1}' i=${i} $infile`
	#<<'COMMENT'
	for cor in `ls ${in_path}/cor_${csta}-???_${tt}.ZZ.sac.${freq}.norm_wha`
	do
		pair=`echo $cor | awk -F "cor_" '{print $2}' | awk -F_ '{print $1}'`
		cor_ZE=`ls ${in_path}/cor_${pair}_${tt}.ZE.sac.${freq}.norm_wha`
		cor_ZN=`ls ${in_path}/cor_${pair}_${tt}.ZN.sac.${freq}.norm_wha`
		rec=`echo $pair | awk -F- '{print$2}'`
		stlo=`awk '$1==rec{print $4}' rec=${rec} UGB21_FALL-Locs_${csta}.txt`
		stla=`awk '$1==rec{print $3}' rec=${rec} UGB21_FALL-Locs_${csta}.txt`
		~/sac_pkg/sacdump $cor | awk '{printf "%.2f %s\n",$1,$2}' > ${pair}_ZZ.tmp
		~/sac_pkg/sacdump $cor_ZE | awk '{printf "%.2f %s\n",$1,$2}' > ${pair}_ZE.tmp
		~/sac_pkg/sacdump $cor_ZN | awk '{printf "%.2f %s\n",$1,$2}' > ${pair}_ZN.tmp
	done
	tmin=-2.0
	tmax=2.0
	
	#=== Start the plot with 0.1 sec inc with amplitude normalization===#
	inc=0.1
	nn=0
	#for ((t=$tmin;t<=$tmax; t+=$inc))
	for t in $(seq $tmin $inc $tmax);
	do
		rm plot.list ppm.list
		ttt=`echo $t | awk '{printf "%.2f",$1}'`
		if [ "$nn" -lt "10" ];then
		nn=`echo 00${nn}`
		elif [ "$nn" -ge "10" ] && [ "$nn" -lt "100" ];then
		nn=`echo 0${nn}`
		fi
		output_ps_file=${ofold}/${nn}_stage_${stage}_${t}_${freq}hz.ps
		for pair in `ls -d CCF/${csta}/${csta}-??? | awk -F/ '{print$3}'`
		do
			rec=`echo $pair | awk -F- '{print$2}'`
			stlo=`awk '$1==rec{print $4}' rec=${rec} UGB21_FALL-Locs_${csta}.txt`
			stla=`awk '$1==rec{print $3}' rec=${rec} UGB21_FALL-Locs_${csta}.txt`
			amp_ZZ=`awk '$1==ttt {print $2}' ttt=${ttt} ${pair}_ZZ.tmp`
			amp_ZE=`awk '$1==ttt {print $2}' ttt=${ttt} ${pair}_ZE.tmp`
			amp_ZN=`awk '$1==ttt {print $2}' ttt=${ttt} ${pair}_ZN.tmp`
			vec=`echo $amp_ZE $amp_ZN | awk '{print sqrt($2**2 + $3**2)}'`
			cos_theta=`echo $amp_ZE $vec | awk '{print ($1)/$2 }'`
			arccos=`echo $amp_ZN $amp_ZE | awk '{print atan2($1,$2)}'`
			if (( $(echo "$amp_ZN >= 0" |bc -l) )); then
				degree=`echo $arccos | awk '{print $1*57.295779513}'`
			else
				degree=`echo $arccos | awk '{print (360 + $1*57.295779513)}'`
			fi
			degree=`echo $arccos | awk '{print $1*57.295779513}'`
			echo $stlo $stla $amp_ZZ >> plot.list
			echo $stlo $stla $degree $vec >> ppm.list
		done

		gmtset ANOT_FONT_SIZE 12
		gmtset BASEMAP_TYPE plain
		gmtset PAPER_MEDIA A3
		gmtset MEASURE_UNIT cm

#		makecpt -Cpolar -T-1.1/1.1/0.001 -Z -Ic > incend.cpt

		cptfile=incend.cpt
		SCA=-JM8i
		REG=-R-110.8301/-110.8292/44.4641/44.4646

		psbasemap $REG $SCA -Gwhite -Ba0.0002f0.0001/a0.0002f0.0001WseN -Xc -Yc -P -K > $output_ps_file
		pscoast $REG $SCA -Df -N2 -K -O -P >> $output_ps_file

		#station amplitude
		awk '{print$1,$2,$3}' plot.list | psxy $SCA $REG -St.7  -W2 -Cincend.cpt -O -K >> $output_ps_file
		#awk '$1==sta{print$4,$3}' sta=${csta} UGB21_FALL-Locs_${csta}.txt | psxy $SCA $REG -Sa.8  -W2 -Gred -O -K >> $output_ps_file
		#Doublet Pool location
		echo -110.829647 44.464345 | psxy $SCA $REG -Sc.8 -W3 -G0/0/0 -O -K >> $output_ps_file
		#awk '{print $1,$2,$3,$4*5}' ppm.list | psxy -J -R -Svt0.005i/0.05i/0.05i -Gblack -W0.5p -K -O >> $output_ps_file
		#awk '{print $1,$2,$3,2*$4}' ppm.list | psxy -J -R -Svt0.005i/0.05i/0.05i -Gblack -W0.5p -K -O >> $output_ps_file

		awk '{print $1,$2,$3,$4}' ppm.list | psxy -J -R -Svt0.005i/0.05i/0.05i -Gblack -W0.5p -K -O >> $output_ps_file


		psscale  -C$cptfile -P -D4i/-0.5/4i/0.5h -B1:'Normalized amp': -K -O  >> $output_ps_file

		pstext -R0/10/0/10 -JX10c -K -O -N -G0 -Y3c -X-1.7c << END >>  $output_ps_file
		5 -1.75 20 0.0 8 CB Minute $stage
END
		pstext -R0/10/0/10 -JX10c -K -O -N -G0 << END >>  $output_ps_file
		5 8 20 0.0 8 CB Second $ttt
END
		convert -trim -density 300 $output_ps_file $output_ps_file.png
		nn=`echo $nn | awk '{print$1+1}'`
		rm $output_ps_file
		#exit
	done
	convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 10 ${ofold}/*_stage_${stage}_*_${freq}hz.ps.png ${ofold}/SNAP_stage_${stage}_${freq}hz.gif
done #end line
#convert -trim -page A4+1+1 -quality 100 -delay 20 $(for i in $(seq -10 1 ); do echo Fig_ppm_single/${eruption}/snr_stage.${i}.ps.png; done) -loop 0 Fig_ppm_single/${eruption}/SNAP_SNR_${eruption}.gif
#convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 10 ${ofold}/${nn}_*_${hz}hz.ps ${ofold}/SNAP_${erup}_${hz}hz.gif
