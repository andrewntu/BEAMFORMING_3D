rm *time_amp_stage*txt
for csta in 562
#for csta in 202
do
	for freq in 1-5
	do
		for stage in $(seq -40 5 50)
		do
			python src/extract_max_value.py $csta $stage $freq OF_SinMei UGB16-Locs-correted_final.txt
		done
	done		
done
