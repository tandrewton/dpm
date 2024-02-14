module load dSQ
#!/bin/bash
numSeeds=2
N_arr=(20 40)
calA0_arr=(1.0)
phi_arr=(0.8)
kl=1.0
ka_arr=(5.0)
kb_arr=(0.01)
#att_arr=(0.0 0.001 0.005 0.01 0.05 0.1)
#att2_arr=(0.0 0.001 0.005 0.01 0.05 0.1)
att_arr=(0.0 0.05 0.1)
att2_arr=(0.0 0.05)
#t_stress_arr=(1.0 10000.0)
t_stress_arr=(10000.0)
v0_arr=(0.1)
gamma_arr=(0)
#gamma_arr=(0 0.001 0.1)
rm joblist_psm.txt
for N in ${N_arr[@]}; do
  for calA0 in ${calA0_arr[@]}; do
    for phi in ${phi_arr[@]}; do
      for ka in ${ka_arr[@]}; do
        for kb in ${kb_arr[@]}; do
          for att in ${att_arr[@]}; do 
            for att2 in ${att2_arr[@]}; do
              for t_stress in ${t_stress_arr[@]}; do
                for v0 in ${v0_arr[@]}; do
                  for gamma in ${gamma_arr[@]}; do
                    echo bash bash/cells/submit_psm.sh 40 30 $calA0 $phi $kl $ka $kb $att $att2 $t_stress $v0 1.0 $gamma 500 pi_ohern,day 0-12:00:00 $numSeeds 1 >> joblist_psm.txt
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

dsq --job-file joblist_psm.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

rm joblist_psm_drawCellSim.txt
for phi in ${phi_arr[@]}; do
  for ka in ${ka_arr[@]}; do
    for kb in ${kb_arr[@]}; do
      for att in ${att_arr[@]}; do 
        for att2 in ${att2_arr[@]}; do
          for t_stress in ${t_stress_arr[@]}; do
            for v0 in ${v0_arr[@]}; do
              for gamma in ${gamma_arr[@]}; do
                echo "module load MATLAB/2023a; matlab -batch \"drawCellSim(\\\"40\\\", \\\"$calA0\\\", \\\"$phi\\\", \\\"$ka\\\", \\\"$kb\\\", \\\"$att\\\", \\\"$att2\\\", \\\"$v0\\\", \\\"$t_stress\\\", \\\"$gamma\\\")\"" >> joblist_psm_drawCellSim.txt
              done
            done
          done
        done
      done
    done
  done
done

**unload matlab before git pull

dsq --job-file joblist_psm_drawCellSim.txt --mem-per-cpu 8g -t 4:00:00 --mail-type NONE --submit --partition pi_ohern,day #--suppress-stats-file -o /dev/null

#dsq --job-file joblist_psm_drawCellSim.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

# send files from mccleary file system to local filesystems
rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

# send files from mccleary output folder (where I store postprocessed files that I ran on the cluster) to local
rsync -rav --inplace --progress --filter '+ */' --filter '+ *sd1.avi' --filter '+ *sd2.avi' --filter '- *.avi' --filter '- *.pos' --filter '- *.tif' --filter '- *' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm/output/ /mnt/c/Users/atata/projects/dpm/output/cells/psm/


rsync -rav --inplace --progress --exclude '*.pos' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /Users/AndrewTon/Documents/YalePhD/projects/dpm/pipeline/cells/. 

rsync -rav --inplace --progress --exclude '*.pos' --exclude '*.tif' --exclude '*.avi' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/home/at965/psm_extracellular_calculation/windowedPhiDataFrame_calA1.0_phi0.8.txt /mnt/c/Users/atata/projects/psm_extracellular_calculation 

salloc -c 4 --mem 16G -t 4:00:00
module load MATLAB/2023a
matlab -nodisplay


# use bash to echo a series of commands that I can copy and paste into a windows terminal to run a python code..
att_arr=(0.001 0.01 0.1)
v0_arr=(0.02 0.04 0.08)
kecm_arr=(0.005 0.05 0.5 5)
for a in ${att_arr[@]}; do
  for v in ${v0_arr[@]}; do
    for e in ${kecm_arr[@]}; do
      echo "python3 minimizedContactsAnalysis.py -a $a -v0 $v -e $e"
    done
  done
done

close all; clear;
calA0_arr = ["1.0"];
att_arr = ["0.001" "0.01" "0.02" "0.05"];
att2_arr = ["0.0" "0.001" "0.01" "0.05"];
phi_arr = ["0.6" "0.8"];
v0_arr = ["0.1"];
gamma_arr = ["0"];
t_stress_arr = ["10000.0"];
ka_arr = ["5.0"];
kb_arr = ["0.01"]; 
%v0_arr = ["0.1"];


for ii=1:length(calA0_arr)
  for jj=1:length(phi_arr)
    for kk=1:length(att_arr)
      for ll=1:length(att2_arr)
        for mm=1:length(v0_arr)
          for nn=1:length(ka_arr)
            for oo=1:length(kb_arr)
              for pp=1:length(t_stress_arr)
                for qq=1:length(gamma_arr)
                  drawCellSim("40", calA0_arr(ii), phi_arr(jj), ka_arr(nn), kb_arr(oo), att_arr(kk), att2_arr(ll), v0_arr(mm), t_stress_arr(pp), gamma_arr(qq))
                end
              end
            end
          end
        end
      end
    end
  end
end


drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.001", "0.1", "10000.0", "0")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.01", "0.1", "1000.0", "0.25")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.01", "0.1", "1000.0", "0.5")

drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0.25")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0.5")