#!/bin/bash
#SBATCH --job-name=parabricks_fq2bam
#SBATCH --partition=gpu
#SBATCH --output=./out/parabricks_fq2bam_%j.out
#SBATCH --error=./out/parabricks_fq2bam_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --nodelist=gpu005,gpu006,gpu007,gpu008,gpu009,gpu010,gpu011,gpu012,gpu013,gpu014,gpu015,gpu016,gpu017,gpu018,gpu019,gpu020,gpu021,gpu022,gpu023,gpu024,gpu025,gpu026,gpu027,gpu028

module load singularity/3.8.3

# Set environment variables for CUDA
export CUDA_VISIBLE_DEVICES=0,1

# Execute each command sequentially

echo "Processing ipsc179"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc179_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc179_2.fq.gz "@RG\\tID:groupipsc179\\tLB:libipsc179\\tPL:ILLUMINA\\tSM:ipsc179\\tPU:groupipsc179" \
    --out-bam ./out/ipsc179_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc196"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc196_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc196_2.fq.gz "@RG\\tID:groupipsc196\\tLB:libipsc196\\tPL:ILLUMINA\\tSM:ipsc196\\tPU:groupipsc196" \
    --out-bam ./out/ipsc196_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc203"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc203_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc203_2.fq.gz "@RG\\tID:groupipsc203\\tLB:libipsc203\\tPL:ILLUMINA\\tSM:ipsc203\\tPU:groupipsc203" \
    --out-bam ./out/ipsc203_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc204"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc204_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc204_2.fq.gz "@RG\\tID:groupipsc204\\tLB:libipsc204\\tPL:ILLUMINA\\tSM:ipsc204\\tPU:groupipsc204" \
    --out-bam ./out/ipsc204_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc206"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc206_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc206_2.fq.gz "@RG\\tID:groupipsc206\\tLB:libipsc206\\tPL:ILLUMINA\\tSM:ipsc206\\tPU:groupipsc206" \
    --out-bam ./out/ipsc206_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc260"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc260_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc260_2.fq.gz "@RG\\tID:groupipsc260\\tLB:libipsc260\\tPL:ILLUMINA\\tSM:ipsc260\\tPU:groupipsc260" \
    --out-bam ./out/ipsc260_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc28"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc28_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc28_2.fq.gz "@RG\\tID:groupipsc28\\tLB:libipsc28\\tPL:ILLUMINA\\tSM:ipsc28\\tPU:groupipsc28" \
    --out-bam ./out/ipsc28_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc399"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc399_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc399_2.fq.gz "@RG\\tID:groupipsc399\\tLB:libipsc399\\tPL:ILLUMINA\\tSM:ipsc399\\tPU:groupipsc399" \
    --out-bam ./out/ipsc399_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc528"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc528_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc528_2.fq.gz "@RG\\tID:groupipsc528\\tLB:libipsc528\\tPL:ILLUMINA\\tSM:ipsc528\\tPU:groupipsc528" \
    --out-bam ./out/ipsc528_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc59"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc59_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc59_2.fq.gz "@RG\\tID:groupipsc59\\tLB:libipsc59\\tPL:ILLUMINA\\tSM:ipsc59\\tPU:groupipsc59" \
    --out-bam ./out/ipsc59_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc608"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc608_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc608_2.fq.gz "@RG\\tID:groupipsc608\\tLB:libipsc608\\tPL:ILLUMINA\\tSM:ipsc608\\tPU:groupipsc608" \
    --out-bam ./out/ipsc608_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc618"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc618_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc618_2.fq.gz "@RG\\tID:groupipsc618\\tLB:libipsc618\\tPL:ILLUMINA\\tSM:ipsc618\\tPU:groupipsc618" \
    --out-bam ./out/ipsc618_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc838"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc838_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc838_2.fq.gz "@RG\\tID:groupipsc838\\tLB:libipsc838\\tPL:ILLUMINA\\tSM:ipsc838\\tPU:groupipsc838" \
    --out-bam ./out/ipsc838_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc864"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc864_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc864_2.fq.gz "@RG\\tID:groupipsc864\\tLB:libipsc864\\tPL:ILLUMINA\\tSM:ipsc864\\tPU:groupipsc864" \
    --out-bam ./out/ipsc864_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc884"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc884_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc884_2.fq.gz "@RG\\tID:groupipsc884\\tLB:libipsc884\\tPL:ILLUMINA\\tSM:ipsc884\\tPU:groupipsc884" \
    --out-bam ./out/ipsc884_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc929"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc929_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc929_2.fq.gz "@RG\\tID:groupipsc929\\tLB:libipsc929\\tPL:ILLUMINA\\tSM:ipsc929\\tPU:groupipsc929" \
    --out-bam ./out/ipsc929_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc92"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc92_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc92_2.fq.gz "@RG\\tID:groupipsc92\\tLB:libipsc92\\tPL:ILLUMINA\\tSM:ipsc92\\tPU:groupipsc92" \
    --out-bam ./out/ipsc92_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

echo "Processing ipsc98"

singularity run --cleanenv --nv \
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
    --in-fq /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc98_1.fq.gz /projectsp/f_ak1833_1/data/X202SC24051878-Z01-F002/01.RawData/ipsc98_2.fq.gz "@RG\\tID:groupipsc98\\tLB:libipsc98\\tPL:ILLUMINA\\tSM:ipsc98\\tPU:groupipsc98" \
    --out-bam ./out/ipsc98_output.bam \
    --bwa-options="-K 10000000" \
    --num-gpus 4 \
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \
    --gpuwrite --gpusort

    