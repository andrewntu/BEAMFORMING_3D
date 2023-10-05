#!/bin/bash

ND=`pwd`
csta=914
freq=1-5hz


for erup in `ls -d ERUPTION_v2/${csta}/${freq}/eruption_*`
#for erup in `ls -d ERUPTION_v2/${csta}/${freq}/eruption_*  | grep -v "eruption_22"`
do
eruption=`echo $erup | awk -F/ '{print$4}'`
for stage in `cat stage_list_1min.txt`
do
#echo "stage " $stage
if [ ! -d "Fig_snr_single" ];then
	mkdir Fig_snr_single/
fi
if [ ! -d "Fig_snr_single/${eruption}" ];then
        mkdir Fig_snr_single/${eruption}
fi
pfile=SNR/$eruption/${csta}_ZZ_SNR_stage.${stage}.${freq}.txt 
Npfile=SNR/$eruption/${csta}_ZN_SNR_stage.${stage}.${freq}.txt 
Epfile=SNR/$eruption/${csta}_ZE_SNR_stage.${stage}.${freq}.txt 
#echo $pfile
if [ "`ls $pfile | wc -l`" -ge "1" ];then
<<'COMMENT'
output_ps_file=Fig_snr_single/${eruption}/snr_stage.$stage.ps
line=`cat $pfile | wc -l`
rm *.tmp
for ((i=1;i<=$line;i++))
do
	pair=`awk 'NR==i{print $1}' i=${i} $pfile | awk -F- '{print$2}'`
	stlo=`awk '$1==pair{print $4}' pair=${pair} UGB21-Locs_sinmei.txt`
	stla=`awk '$1==pair{print $3}' pair=${pair} UGB21-Locs_sinmei.txt`
	snr=`awk 'NR==i{print $6}' i=${i} $pfile`
	Nsnr=`awk 'NR==i{print $6}' i=${i} $Npfile`
	Esnr=`awk 'NR==i{print $6}' i=${i} $Epfile`
#	echo ${pair}.ZZ.stage.${stage}.sac.norm $stlo $stla
	echo $stlo $stla $snr >> plot.tmp
	echo $stlo $stla $Nsnr >> Nplot.tmp
	echo $stlo $stla $Esnr >> Eplot.tmp
done
#<<'COMMENT'
gmtset ANOT_FONT_SIZE 10
gmtset BASEMAP_TYPE plain
gmtset PAPER_MEDIA A3
gmtset MEASURE_UNIT cm

makecpt -Cpanoply -T0/8/0.01 -Z -Ic > incend.cpt

cptfile=incend.cpt
SCA=-JM5i
REG=-R-110.831/-110.8264/44.4593/44.462

psbasemap $REG $SCA -Gwhite -B0.001/0.001WsNe -Xc -Y12i -P -K > $output_ps_file

#vertical
cat plot.tmp | grep -v "nan" | psxy $SCA $REG -Sc.35  -C$cptfile -O -K >> $output_ps_file
#Old Faithful location
#echo -110.828211897200987 44.460437153479248 | psxy $SCA $REG -Sa.8 -W3 -G0/0/0 -O -K >> $output_ps_file


psbasemap $REG $SCA -Gwhite -B0.001/0.001WsNe -Y-5i -P -K -O >> $output_ps_file
#vertical
cat Nplot.tmp | grep -v "nan" | psxy $SCA $REG -Sc.35  -C$cptfile -O -K >> $output_ps_file
#Old Faithful location
#echo -110.828211897200987 44.460437153479248 | psxy $SCA $REG -Sa.8 -W3 -G0/0/0 -O -K >> $output_ps_file


psbasemap $REG $SCA -Gwhite -B0.001/0.001WsNe -Y-5i -P -K -O >> $output_ps_file
#vertical
cat Eplot.tmp | grep -v "nan" | psxy $SCA $REG -Sc.35  -C$cptfile -O -K >> $output_ps_file
#Old Faithful location
#echo -110.828211897200987 44.460437153479248 | psxy $SCA $REG -Sa.8 -W3 -G0/0/0 -O -K >> $output_ps_file



psscale  -C$cptfile -P -D6/-0.5/10/0.5h -B1:'SNR': -K -O  >> $output_ps_file

pstext -R0/10/0/10 -JX10c -K -O -N -G0 -Y3c -X-1.7c << END >>  $output_ps_file
3.45 -1.75 12 0.0 7 CB  Minute $stage
END
#exit
#COMMENT
convert -trim -density 300 $output_ps_file $output_ps_file.png
rm $output_ps_file
COMMENT
fi
done
convert -trim -page A4+1+1 -quality 100 -delay 20 $(for i in $(seq -40 1 70); do echo Fig_snr_single/${eruption}/snr_stage.${i}.ps.png; done) -loop 0 Fig_snr_single/${eruption}/SNAP_SNR_${eruption}.gif
done #stage
