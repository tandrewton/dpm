module load dSQ
#!/bin/bash
numSeeds=1
calA0=(1.0)
phi_arr=(0.8)
kl=1.0
ka_arr=(5.0)
kb_arr=(0.1)
att_arr=(0.005 0.05)
att2_arr=(0.05)
v0_arr=(0.0 0.05 0.1)
gamma_arr=(0 0.25 0.5)
rm joblist_psm_att_v0.txt
for phi in ${phi_arr[@]}; do
  for ka in ${ka_arr[@]}; do
    for kb in ${kb_arr[@]}; do
      for att in ${att_arr[@]}; do 
        for att2 in ${att2_arr[@]}; do
          for v0 in ${v0_arr[@]}; do
            for gamma in ${gamma_arr[@]}; do
              echo bash bash/cells/submit_psm.sh 20 30 $calA0 $phi $kl $ka $kb $att $att2 0 $v0 1.0 $gamma 100 pi_ohern,day 0-12:00:00 $numSeeds 1 >> joblist_psm_att_v0.txt
            done
          done
        done
      done
    done
  done
done

dsq --job-file joblist_psm_att_v0.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

bash bash/cells/submit_psm.sh 40 20 1.05 0.9 0 0.0 0.05 50.0 1.0 10.0 1000 pi_ohern,day 0-4:00:00 1 1

rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

rsync -rav --inplace --progress at965@transfer-mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /Users/AndrewTon/Documents/YalePhD/projects/dpm/pipeline/cells/. 

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
%att_arr = ["0.005" "0.05"];
att2_arr = ["0.05"];
phi_arr = ["0.8"];
%ka_arr = ["1.0" "5.0"];
%kb_arr = ["0.01" "0.1"];
%v0_arr = ["0.0" "0.05" "0.1"];
%gamma_arr = ["0" "0.25" "0.5" "1.0"];

att_arr = ["0.05"];
ka_arr = ["5.0"];
kb_arr = ["0.1"];
v0_arr = ["0.1"];
gamma_arr = ["0"];


for ii=1:length(calA0_arr)
  for jj=1:length(phi_arr)
    for kk=1:length(att_arr)
      for ll=1:length(att2_arr)
        for mm=1:length(v0_arr)
          for nn=1:length(ka_arr)
            for oo=1:length(kb_arr)
              for pp=1:length(gamma_arr)
                drawCellSim("20", calA0_arr(ii), phi_arr(jj), ka_arr(nn), kb_arr(oo), att_arr(kk), att2_arr(ll), v0_arr(mm), gamma_arr(pp))
              end
            end
          end
        end
      end
    end
  end
end