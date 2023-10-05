#!/bin/bash

ND=`pwd`
if [ $# != 3 ]; then
  echo "Usage format.csh: [src] [tt] [hz]"
  echo
  exit
fi

src=$1
tt=$2
hz=$3
#<<'COMMENT'
whf=amp_${tt}_${hz}hz
sta_file=UGB21-Locs_sinmei_final.txt
sampling_tt=`echo $hz | awk -F"-" '{print$2}' | awk '{print 1/$1}'`
#sampling_tt=0.1
#<<'COMMENT'
rm -rf $whf $whf.sort 
#<<'COMMENT'
for wave in `ls Data_CCF_Nan/${src}/${hz}hz/${src}-???.ZZ.stage.${tt}.sac.norm`
do
name=`echo $wave | awk -F/ '{print$4}' | awk -F. '{print$1}'`
ZZ_name=`echo $wave | awk -F/ '{print$4}'`
ZE_wave=`ls Data_CCF_Nan/${src}/${hz}hz/${name}.ZE.stage.${tt}.sac.norm`
ZN_wave=`ls Data_CCF_Nan/${src}/${hz}hz/${name}.ZN.stage.${tt}.sac.norm`
#echo $name $ZE_wave $ZN_wave
#exit
ZE_name=`echo $ZE_wave | awk -F/ '{print$4}'`
ZN_name=`echo $ZN_wave | awk -F/ '{print$4}'`
rm Figure/sac_ascii/${ZZ_name}_${hz}.txt Figure/sac_ascii/${ZN_name}_${hz}.txt Figure/sac_ascii/${ZE_name}_${hz}.txt
saclst depmax f $wave | awk '{print $2}' >> $whf
~/sac_pkg/sacdump $wave > Figure/sac_ascii/${ZZ_name}_${hz}.txt
~/sac_pkg/sacdump $ZE_wave > Figure/sac_ascii/${ZE_name}_${hz}.txt
~/sac_pkg/sacdump $ZN_wave > Figure/sac_ascii/${ZN_name}_${hz}.txt
done
#COMMENT
cat $whf | sort -k1 -g -r > $whf.sort
#COMMENT


gmtset ANOT_FONT_SIZE 10
gmtset BASEMAP_TYPE plain
gmtset PAPER_MEDIA A3
gmtset MEASURE_UNIT cm


SCA=-JM10i
#REG=-R-110.831/-110.8264/44.4593/44.462
REG=-R-110.829/-110.827/44.46/44.461
makecpt -Cjet -T-1.5/1.5/0.001 -Z > tt.cpt
#makecpt -Cjet -T-1.0/1.0/0.001 > amp.cpt
makecpt -Cpolar -T-1.0/1.0/0.001 > amp.cpt
#makecpt -Cjet -T-1.0/1.0/0.00001 -Z > amp.cpt
nn=0
max_amp=`head -n 1 amp_${tt}_${hz}hz.sort`

if [ ! -d "Figure/${src}" ];then
        mkdir Figure/${src}
fi

if [ ! -d "Figure/${src}/erup_${tt}" ];then
        mkdir Figure/${src}/erup_${tt}
fi
rm Figure/${src}/erup_${tt}/*erup_${tt}_${hz}hz*jpg
#for i in $(seq -2.0 $sampling_tt 2.0)
#for i in $(seq -1.0 $sampling_tt 1.0)
#for i in $(seq -0.5 $sampling_tt 0.5)
for i in $(seq -0.5 0.01 0.5)
#for i in $(seq -0.0 0.01 0.01)
do
nn=`echo $nn | awk '{print$1+1}'`
if [ "$nn" -lt "10" ];then
nn=`echo 00${nn}`
elif [ "$nn" -ge "10" ] && [ "$nn" -lt "100" ];then
nn=`echo 0${nn}`
fi
rm plot.list
rm ppm.list
output_ps_file=Figure/${src}/erup_${tt}/${nn}_${i}_erup_${tt}_${hz}hz.ps
echo "writing $output_ps_file file"
if [ -f $output_ps_file ];then
rm $output_ps_file
fi
#gmtBASEMAP_TYPE fancy
psbasemap $REG $SCA -Gwhite -Ba0.001f0.001/a0.001f0.001WSen:". ${i} sec": -Xc -Yc -P -V -K > $output_ps_file
#psbasemap $REG $SCA -Gwhite -Ba0.001:"TT ${i} sec":/f0.05WSen -Xc -Yc -P -V -K > $output_ps_file
for file in `ls Figure/sac_ascii/${src}-???.ZZ.stage.${tt}.sac.norm_${hz}.txt`
do

file_name=`echo $file | awk -F/ '{print$3}' | awk -F. '{print$1}'`
ZE_file=Figure/sac_ascii/${file_name}.ZE.stage.${tt}.sac.norm_${hz}.txt
ZN_file=Figure/sac_ascii/${file_name}.ZN.stage.${tt}.sac.norm_${hz}.txt

sta=`echo $file | awk -F/ '{print$3}' | awk -F"-" '{print$2}' | awk -F. '{print$1}'`
#amp=`awk '$1==i{print $2/max_amp}' i=${i} max_amp=${max_amp} $file`
#amp=`awk '{printf "%.2f %s\n",$1,$2}' $file | awk '$1==i{print $2/max_amp}' i=${i} max_amp=${max_amp}`
amp=`awk '{printf "%.2f %s\n",$1,$2}' $file | awk '$1==i{print $2}' i=${i}`
amp_ZE=`awk '{printf "%.2f %s\n",$1,$2}' $ZE_file | awk '$1==i{print $2}' i=${i}`
amp_ZN=`awk '{printf "%.2f %s\n",$1,$2}' $ZN_file | awk '$1==i{print $2}' i=${i}`

vec=`echo $amp_ZE $amp_ZN | awk '{print sqrt($1**2 + $2**2)}'`
# (ZE,ZN)(1,0) = cos
cos_theta=`echo $amp_ZE $vec | awk '{print ($1)/$2 }'`
#arccos=`echo $cos_theta | awk '{print atan2()}'`
arccos=`echo $amp_ZN $amp_ZE | awk '{print atan2($1,$2)}'`
#echo $arccos
#<<'COMMENT'
if (( $(echo "$amp_ZN >= 0" |bc -l) )); then
	degree=`echo $arccos | awk '{print $1*57.295779513}'`
else
	degree=`echo $arccos | awk '{print (360 + $1*57.295779513)}'`
fi
#COMMENT
degree=`echo $arccos | awk '{print $1*57.295779513}'`
#echo $amp_ZN $amp_ZE $degree
stlo=`awk '$1==sta{print$4}' sta=${sta} $sta_file`
stla=`awk '$1==sta{print$3}' sta=${sta} $sta_file`
#echo $stlo $stla $amp
echo $stlo $stla $amp >> plot.list
echo $stlo $stla $degree $vec >> ppm.list 
#echo $stlo $stla $amp_ZE $amp_ZN >> ppm.list 
done
awk '{print$1,$2,$3}' plot.list | psxy $SCA $REG -St.5  -W2 -Camp.cpt -O -K >> $output_ps_file
echo -110.828211897200987 44.460437153479248 | psxy $SCA $REG -Sc.8 -W3 -G0/0/0 -O -K >> $output_ps_file
awk '$1==sta{print$4,$3}' sta=${src} $sta_file | psxy $SCA $REG -Sa.8  -W2 -Gred -O -K >> $output_ps_file
awk '{print $1,$2,$3,$4*5}' ppm.list | psxy -J -R -Svt0.005i/0.05i/0.05i -Gblack -W0.5p -K -O >> $output_ps_file
#awk '{print $1,$2,$3*0.75,$4*0.75}' ppm.list | psvelo -J -R −H2 -W0.5p -Gblack  -K -O >> $output_ps_file

psscale -Camp.cpt -P -D12.5/-1/15/1/h -Ba0.2:"Amplitude":/f0.05 -O -K  >> $output_ps_file
psxy -R -J -O -T >> $output_ps_file
convert -density 300 -trim $output_ps_file $output_ps_file'.jpg'
rm $output_ps_file
#exit
done
rm *.cpt
#COMMENT
if [ "${hz}" == "1-5" ]; then
#convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 100 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/SNAP_erup_${tt}_${hz}hz.gif
convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 10 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/${src}/SNAP_erup_${tt}_${hz}hz.gif
#convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 20 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/SNAP_erup_${tt}_${hz}hz_d20.gif
else
#convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 20 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/SNAP_erup_${tt}_${hz}hz.gif
convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 10 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/${src}/SNAP_erup_${tt}_${hz}hz.gif
#convert -trim -page A4+1+1 -quality 100  -loop 0 -delay 20 Figure/${src}/erup_${tt}/*_erup_${tt}_${hz}hz.ps.jpg Figure/SNAP_erup_${tt}_${hz}hz_d20.gif
fi
