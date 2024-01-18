module load dSQ
#!/bin/bash
numSeeds=10
calA0=(1.0)
phi_arr=(0.8)
kl=1.0
ka_arr=(5.0)
kb_arr=(0.1)
att_arr=(0.001 0.05)
att2_arr=(0.001 0.05)
t_stress_arr=(1.0 10000.0)
v0_arr=(0.0 0.1)
gamma_arr=(0 0.25 0.5)
rm joblist_psm_att_v0.txt
for phi in ${phi_arr[@]}; do
  for ka in ${ka_arr[@]}; do
    for kb in ${kb_arr[@]}; do
      for att in ${att_arr[@]}; do 
        for att2 in ${att2_arr[@]}; do
          for t_stress in ${t_stress_arr[@]}; do
            for v0 in ${v0_arr[@]}; do
              for gamma in ${gamma_arr[@]}; do
                echo bash bash/cells/submit_psm.sh 40 30 $calA0 $phi $kl $ka $kb $att $att2 $t_stress $v0 1.0 $gamma 200 pi_ohern,day 0-12:00:00 $numSeeds 1 >> joblist_psm_att_v0.txt
              done
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
att_arr = ["0.001" "0.05"];
att2_arr = ["0.001" "0.05"];
phi_arr = ["0.8"];
v0_arr = ["0.0" "0.1"];
gamma_arr = ["0" "0.25" "0.5"];
t_stress_arr = ["1.0" "10000.0" ];
ka_arr = ["5.0"];
kb_arr = ["0.1"];
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

# can I make movies on the cluster?
# turn this into a bash script and use dsq to submit on cluster. 

#!/bin/bash
#SBATCH --job-name myjob
#SBATCH --cpus-per-task 4
#SBATCH --mem 18G
#SBATCH -t 8:00:00

module load MATLAB/2023a
# assuming you have your_script.m in the current directory
matlab -batch "drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.001", "0.1", "10000.0", "0")"

# if using MATLAB older than R2019a
# matlab -nojvm -nodisplay -nosplash < your_script.m


>> profiles = VideoWriter.getProfiles()
  Summary of installed VideoWriter profiles:

          Name                                     Description
    ---------------- -----------------------------------------------------------------------
    Archival         Video file compression with JPEG 2000 codec with lossless mode enabled.
    Grayscale AVI    An AVI file with Grayscale Video Data
    Indexed AVI      An AVI file with Indexed Video Data
    Motion JPEG 2000 Video file compression with JPEG 2000 codec.
    Motion JPEG AVI  An AVI file with Motion JPEG compression
    Uncompressed AVI An AVI file with uncompressed RGB24 video data

Name=psm_calA01.0_phi0.8_tm1000.0_v00.1_t_abp1.0_gamma0.25_kl1.0_ka5.0_kb0.1_N40_dur100_att0.001_att20.01_start1_end1 
Name=psm_calA01.0_phi0.8_tm1000.0_v00.1_t_abp1.0_gamma0_kl1.0_ka5.0_kb0.1_N40_dur100_att0.05_att20.001_start1_end1

drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.001", "0.1", "10000.0", "0")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.01", "0.1", "1000.0", "0.25")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.01", "0.1", "1000.0", "0.5")

drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0.25")
drawCellSim("40", "1.0", "0.8", "5.0", "0.1", "0.001", "0.05", "0.0", "1.0", "0.5")