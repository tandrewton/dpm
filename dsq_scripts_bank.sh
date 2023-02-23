module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in ${att_arr[@]}; do
    for om in 0.01; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_att.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_att.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in 0.1; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_om.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_om.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in ${t_stress_arr[@]}; do
  for att in 0.1; do
    for om in 0.01; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_tau_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_tau_ka.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_att_om.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_om.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in 0.01; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_att_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_ka.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null

module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in 100000.0; do
  for att in 0.1; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_om_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_om_ka.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null


module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in 100000.0; do
  for att in 0.1; do
    for om in 0.01; do
      for kl in ${kl_arr[@]}; do
        for ka in ${ka_arr[@]}; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_kl_ka.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_kl_ka.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null

# runs for presentation
module load dSQ
# testing stress relaxation, attraction, boundaries on
#bd0 P
#!/bin/bash
numSeeds=2
t_stress_arr=(1.0 5.0 25.0 125.0 100000.0)
att_arr=(0.05 0.1 0.15 0.2 0.25)
om_arr=(0.001 0.005 0.01 0.05)
kl_arr=(0.1 0.5 1.0 5.0 10.0)
ka_arr=(0.1 0.5 1.0 5.0 10.0)
rm joblist_PS.txt
for t_stress in 100000.0; do
  for att in ${att_arr[@]}; do
    for om in ${om_arr[@]}; do
      for kl in 1.0; do
        for ka in 1.0; do
            echo bash bash/epi2D/submit_laserAblation.sh 50 30 3 1.10 0.94 0.85 $kl $ka 0.01 $att $om 4.0 1.0 4.0 1.0 0.0 $t_stress 0 1 1000 pi_ohern,day,scavenge 0-24:00:00 $numSeeds 1 >> joblist_PS_att_om.txt
        done
      done
    done
  done
done

dsq --job-file joblist_PS_att_om.txt --mem-per-cpu 4g -t 24:00:00 --mail-type NONE --submit --suppress-stats-file  -o /dev/null

