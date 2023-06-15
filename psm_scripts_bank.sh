
module load dSQ
#!/bin/bash
numSeeds=1
att_arr=(0 0.01 0.05 0.1)
v0_arr=(0.01 0.05 0.1)
koff_arr=(0.1 1.0 10.0 100.0 1000.0)
rm joblist_PS_tau_ka.txt
for att in ${att_arr[@]}; do 
  for v0 in ${v0_arr[@]}; do
    for koff in ${koff_arr[@]}; do
        echo bash bash/epi2D/submit_psm.sh 40 20 1.05 0.9 $att 0.0 $v0 50.0 $koff 1000 day 0-24:00:00 $numSeeds 1 >> joblist_psm_att_v0_koff.txt
    done
  done
done

dsq --job-file joblist_psm_att_v0_koff.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null