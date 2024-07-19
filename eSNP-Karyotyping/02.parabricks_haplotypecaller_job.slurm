#!/bin/bash
#SBATCH --job-name=parabricks_haplotypecaller
#SBATCH --partition=gpu
#SBATCH --output=./out/parabricks_haplotypecaller_%j.out
#SBATCH --error=./out/parabricks_haplotypecaller_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --mem=60G
#SBATCH --time=24:00:00

module load singularity/3.8.3

# Set environment variables for CUDA
export CUDA_VISIBLE_DEVICES=0,1,2,3

# Execute each command sequentially

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc179.bam \
    --out-variants ./out/ipsc179.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc196.bam \
    --out-variants ./out/ipsc196.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc203.bam \
    --out-variants ./out/ipsc203.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc204.bam \
    --out-variants ./out/ipsc204.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc206.bam \
    --out-variants ./out/ipsc206.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc260.bam \
    --out-variants ./out/ipsc260.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc28.bam \
    --out-variants ./out/ipsc28.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc28_output.bam \
    --out-variants ./out/ipsc28_output.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc399.bam \
    --out-variants ./out/ipsc399.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc528.bam \
    --out-variants ./out/ipsc528.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc59.bam \
    --out-variants ./out/ipsc59.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc608.bam \
    --out-variants ./out/ipsc608.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc618.bam \
    --out-variants ./out/ipsc618.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc838.bam \
    --out-variants ./out/ipsc838.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc864.bam \
    --out-variants ./out/ipsc864.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc884.bam \
    --out-variants ./out/ipsc884.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc92.bam \
    --out-variants ./out/ipsc92.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc929.bam \
    --out-variants ./out/ipsc929.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-bam ./out/ipsc98.bam \
    --out-variants ./out/ipsc98.vcf \
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10

    