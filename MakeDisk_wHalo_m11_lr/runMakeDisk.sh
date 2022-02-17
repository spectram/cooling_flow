#!/bin/bash

#SBATCH --job-name=makeDisk
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=4:00:00
#SBATCH --output=/mnt/home/jstern/cooling_flow/MakeDisk_wHalo_m11_lr/job_output/%j.out
#SBATCH --error=/mnt/home/jstern/cooling_flow/MakeDisk_wHalo_m11_lr/job_output/%j.err
#SBATCH --mail-user=sternnaty@gmail.com
#SBATCH --mail-type=fail
#SBATCH --mail-type=end

./MakeHubbleType ~/ceph/ICs/m11_no_halo_150_res2e3_fgas02.ic


