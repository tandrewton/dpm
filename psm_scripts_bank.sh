module load dSQ
#!/bin/bash
numSeeds=1
calA0=(1.0)
phi_arr=(0.74 0.63 0.53)
att_arr=(0 0.001 0.01 0.1)
v0_arr=(0.005)
kecm_arr=(0.006 0.06 0.6)
koff_arr=(1.0 10.0 100.0)
rm joblist_psm_att_v0_koff.txt
for phi in ${phi_arr[@]}; do
  for att in ${att_arr[@]}; do 
    for v0 in ${v0_arr[@]}; do
      for koff in ${koff_arr[@]}; do
          for kecm in ${kecm_arr[@]}; do
              echo bash bash/cells/submit_psm.sh 40 20 $calA0 $phi $att 0.0 $v0 1.0 $kecm $koff 400 pi_ohern,day 0-4:00:00 $numSeeds 1 >> joblist_psm_att_v0_koff.txt
          done
      done
    done
  done
done

dsq --job-file joblist_psm_att_v0_koff.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

bash bash/cells/submit_psm.sh 40 20 1.05 0.9 0 0.0 0.05 50.0 1.0 10.0 1000 pi_ohern,day 0-4:00:00 1 1

rsync -rav --inplace --progress at965@transfer-grace.hpc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

close all; clear;
calA0_arr = ["1.0"];
att_arr = ["0" "0.001" "0.01" "0.1"];
phi_arr = ["0.74" "0.63" "0.53"]
v0_arr = ["0.005"];
k_ecm_arr = ["0.006" "0.06" "0.6"];
k_off_arr = ["1.0" "10.0" "100.0"];
for ii=1:length(calA0_arr)
    for jj=1:length(phi_arr)
        for kk=1:length(att_arr)
            for ll=1:length(v0_arr)
                for mm=1:length(k_ecm_arr)
                  for nn=1:length(k_off_arr)
                    drawCellSim("40", calA0_arr(ii), phi_arr(jj), att_arr(kk), v0_arr(ll), k_ecm_arr(mm), k_off_arr(nn))
                  end
                end
            end
        end
    end
end