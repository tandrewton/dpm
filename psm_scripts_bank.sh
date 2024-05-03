module load dSQ
#!/bin/bash
numSeeds=10
N_arr=(40)
calA0_arr=(1.0)
phi_arr=(0.8)
kl_arr=(0.2)
ka=2.5
kb_arr=(0.01)
#att_arr=(0.005 0.01 0.03 0.05)
#att2_arr=(0 0.0005 0.001 0.005 0.01 0.05 0.1)
att_arr=(0.01 0.03)
att2_arr=(0 0.005 0.01 0.05)
#kecm_arr=(0.001 0.005 0.01 0.05)
t_stress_arr=(10000.0)
v0_arr=(0.1)
gamma_arr=(0)
kon_arr=(1.0)
koff_arr=(0.1 0.5 1.0 10.0 100.0)
calcMinPos=1
rm joblist_psm.txt
for N in ${N_arr[@]}; do
  for calA0 in ${calA0_arr[@]}; do
    for phi in ${phi_arr[@]}; do
      for kl in ${kl_arr[@]}; do
        for kb in ${kb_arr[@]}; do
          for att in ${att_arr[@]}; do 
            for att2 in ${att2_arr[@]}; do
              for t_stress in ${t_stress_arr[@]}; do
                for v0 in ${v0_arr[@]}; do
                  for gamma in ${gamma_arr[@]}; do
                    for k_on in ${kon_arr[@]}; do
                      for k_off in ${koff_arr[@]}; do
                        k_ecm=$att2
                        #for k_ecm in ${kecm_arr[@]}; do
                        echo bash bash/cells/submit_psm.sh $N 30 $calA0 $phi $kl $ka $kb $att $att2 $t_stress $v0 1.0 $gamma $k_on $k_off $k_ecm $calcMinPos 300 pi_ohern,day 0-24:00:00 $numSeeds 1 >> joblist_psm.txt
                        #done
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
  done
done

dsq --job-file joblist_psm.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

rm joblist_psm_drawCellSim.txt
for N in ${N_arr[@]}; do
  for calA0 in ${calA0_arr[@]}; do
    for phi in ${phi_arr[@]}; do
      for kl in ${kl_arr[@]}; do
        for kb in ${kb_arr[@]}; do
          for att in ${att_arr[@]}; do 
            for att2 in ${att2_arr[@]}; do
              for t_stress in ${t_stress_arr[@]}; do
                for v0 in ${v0_arr[@]}; do
                  for gamma in ${gamma_arr[@]}; do
                    for k_on in ${kon_arr[@]}; do
                      for k_off in ${koff_arr[@]}; do
                        k_ecm=$att2
                        #for k_ecm in ${kecm_arr[@]}; do
                        echo "module load MATLAB/2023a; matlab -batch \"drawCellSim(\\\"$N\\\", \\\"$calA0\\\", \\\"$phi\\\", \\\"$kl\\\", \\\"$kb\\\", \\\"$att\\\", \\\"$att2\\\", \\\"$v0\\\", \\\"$t_stress\\\", \\\"$gamma\\\", \\\"$k_on\\\", \\\"$k_off\\\", \\\"$k_ecm\\\",$numSeeds)\"" >> joblist_psm_drawCellSim.txt
                        #done
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
  done
done

dsq --job-file joblist_psm_drawCellSim.txt --mem-per-cpu 8g -t 8:00:00 --mail-type NONE --submit --partition pi_ohern,day #--suppress-stats-file -o /dev/null

cd ~/psm_extracellular_calculation/
rm windowedPhiDataFrame_calA*
salloc -c 4 --mem 16G -t 4:00:00
module load MATLAB/2023a
matlab -nodisplay

delete windowedPhiDataFrame.txt
copy paste the contents of windowedFractionalPhiDP
then run cellShapeAnalysis.py


#dsq --job-file joblist_psm_drawCellSim.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

# send files from mccleary file system to local filesystems
#rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

# send files from mccleary output folder (where I store postprocessed files that I ran on the cluster) to local
# expecting speed.csv, shape.csv, and sd1.avi files
rsync -rav --inplace --progress --filter '+ */' --filter '+ *sd1.avi' --filter '- *.avi' --filter '+ *sd1.mp4' --filter '- *.mp4' --filter '+ *.csv' --filter '- *.pos'  --filter '+ *fr60.tif' --filter '- *.tif' --filter '+ */' --filter '- *' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm/output/ /mnt/c/Users/atata/projects/dpm/output/cells/psm/
rsync -rav --inplace --progress --filter '+ */' --filter '+ *sd1.avi' --filter '- *.avi' --filter '+ *.csv' --filter '- *.pos'  --filter '+ *fr95.tif' --filter '- *.tif' --filter '+ */' --filter '- *' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm/output/ /Users/AndrewTon/Documents/YalePhD/projects/dpm/output/cells/psm/

# expecting .xstream, .xminstream, .shapestream, etc
rsync -rav --inplace --progress --exclude '*.pos' --exclude '*.tif' --exclude '*.avi' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 
rsync -rav --inplace --progress --exclude '*.pos' --exclude '*.tif' --exclude '*.avi' at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /Users/AndrewTon/Documents/YalePhD/projects/dpm/pipeline/cells/. 

# expecting windowedPhiDataFrame*.txt
rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/home/at965/psm_extracellular_calculation/windowedPhiDataFrame_calA1.0_phi0.8.txt /mnt/c/Users/atata/projects/psm_extracellular_calculation 
rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/home/at965/psm_extracellular_calculation/windowedPhiDataFrame_calA1.0_phi0.8.txt /Users/AndrewTon/Documents/YalePhD/projects/ZebrafishSegmentation/psm_extracellular_calculation 


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
phi_arr = ["0.8"];
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

for i in *.avi; do ffmpeg -i "$i" -c:v libx264 -preset slow  -profile:v high -level:v 4.0 -pix_fmt yuv420p -crf 22 -codec:a aac "$(basename "$i" .avi)".mp4  ; done
ffmpeg -i input.avi -c:v copy -c:a copy OUTPUT.mp4