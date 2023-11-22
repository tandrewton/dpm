module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_tau_att.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in ${att_arr[@]}; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_att.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_att.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_tau_om.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in 0.1; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_om.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_om.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
# t_stress_arr=(1.0 2.0 4.0 8.0 16.0 32.0 64.0 128.0 256.0 512.0 1024.0 100000.0)
t_stress_arr=(2.0 8.0 32.0 128.0 512.0 100000.0)
#t_stress_arr=(1.0 4.0 16.0 64.0 256.0 1024.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.2 0.6 1.0 1.4 1.8 2.2 2.6 3.0 3.4 3.8 4.2 10.0) 
rm joblist_PS_tau_ka.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in 0.1; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1400 day 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_ka.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_att_om.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_att_om.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_om.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_att_ka.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_att_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_ka.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_att_ka.txt
for t_stress in 100000.0; do
  for att in 0.1; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_att_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_ka.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null
 
module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_kl_ka.txt
for t_stress in 100000.0; do
  for att in 0.1; do
    for om in 1.0; do
      for kl in ${kl_arr[@]}; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_kl_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_kl_ka.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
rm joblist_PS_tau_att.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in ${att_arr[@]}; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka 0.001 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_att.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_att.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
kb_arr=(0.001 0.01 0.1)
rm joblist_PS_tau_kb.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in 0.1; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in 1.0; do
          for kb in ${kb_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka $kb $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_kb.txt
          done
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_kb.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=10
t_stress_arr=(1.0 5.0 25.0 125.0 250.0 625.0 100000.0)
att_arr=(0.01 0.02 0.05 0.1 0.2)
om_arr=(0.1 0.5 1.0 5.0)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 2.5 5.0 7.5 10.0)
kb_arr=(0.001 0.01 0.1)
rm joblist_PS_att_kb.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in 1.0; do
      for kl in 1.0; do
        for ka in 1.0; do
          for kb in ${kb_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.20 0.94 0.85 $kl $ka $kb $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 day 0-24:00:00 $numSeeds 1 >> joblist_PS_att_kb.txt
          done
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_kb.txt --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

# for wing disc and embryo production simulations
module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=50
#t_stress_arr=(4.8 9.6 19.2 76.8 307.2 1228.8 4915.2 9830.4 39321.6) #full range
#t_stress_arr=(19.2 38.4 76.8) #embryo
#t_stress_arr=(4915.2 9830.4) # wing disc
t_stress_arr=(307.2 614.4)

#ka_arr=(0.5 1.0 2.0 4.0 8.0 16.0 32.0 40.0)
ka_arr=(32.0)

#tauRatio_arr=(0.1 0.2 0.3 0.4)
#tauRatio_arr=(0.1 0.4)
tauRatio_arr=(0 0.1 0.2 0.3 0.4)
k_ps=(4.0)
rm joblist_PS_tau_ka.txt
for t_stress in ${t_stress_arr[@]}; do 
  for taur in ${tauRatio_arr[@]}; do 
    for att in 0.1; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            if [ "$(echo "$ka < 23.0" | bc)" -eq 1 ]; then 
              clusters="pi_ohern,day"
              runtime="0-24:00:00"
            else
              clusters="pi_ohern,week"
              runtime="2-24:00:00"
            fi
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 5 1.20 0.94 0.99 $kl $ka 0.01 $att 1.0 4.0 $k_ps 4.0 1.0 0.0 $t_stress $taur 1 2000 $clusters $runtime $numSeeds 1 >> joblist_PS_tau_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_ka.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

./main/epi2D/laserAblation.o 30 20  5 1.20 0.94 0.85 1.0 16.0 0.01 0.1 1.0  4.0  4.0  4.0 1.0  0.0   9830.4  0  1  1 500  test