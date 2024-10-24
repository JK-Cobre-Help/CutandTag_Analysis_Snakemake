# CutandTag_Analysis_Snakemake

## 1. Clone repository
```
git clone https://github.com/JK-Cobre-Help/CutandTag_Analysis_Snakemake.git
```
## 2. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
## 3. Modify Samples file
```
vim samples.csv
```
## 4. Dry Run
```
snakemake -npr
```
## 5. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 30 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```
