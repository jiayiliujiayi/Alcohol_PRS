import os
import argparse

def generate_slurm_script(input_dir, output_dir, output_filename, read_group_prefix):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Find all .fq.gz files in the input directory
    fq_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.fq.gz')])

    # Assuming paired-end reads
    pairs = {}
    for file in fq_files:
        prefix = '_'.join(file.split('_')[:-1])  # Modify this if file naming convention changes
        if prefix not in pairs:
            pairs[prefix] = []
        pairs[prefix].append(file)

    # Check for any unpaired files
    for key, files in pairs.items():
        if len(files) != 2:
            raise ValueError(f"Unpaired files found for prefix {key}: {files}")

    # Generate commands for each pair
    commands = []
    for index, (prefix, files) in enumerate(pairs.items()):
        bam_path = os.path.join(output_dir, f"{prefix}_output.bam")
        command = f"""
echo "Processing {prefix}"

singularity run --cleanenv --nv \\
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \\
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun fq2bamfast \\
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \\
    --in-fq {os.path.join(input_dir, files[0])} {os.path.join(input_dir, files[1])} "@RG\\\\tID:{read_group_prefix}{prefix}\\\\tLB:lib{prefix}\\\\tPL:ILLUMINA\\\\tSM:{prefix}\\\\tPU:{read_group_prefix}{prefix}" \\
    --out-bam {bam_path} \\
    --bwa-options="-K 10000000" \\
    --num-gpus 4 \\
    --bwa-nstreams 2 --bwa-cpu-thread-pool 3  \\
    --gpuwrite --gpusort
"""
        commands.append(command)

    # Create the SLURM script content
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=parabricks_fq2bam
#SBATCH --partition=gpu
#SBATCH --output={output_dir}/parabricks_fq2bam_%j.out
#SBATCH --error={output_dir}/parabricks_fq2bam_%j.err
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
{''.join(commands)}
    """

    # Write the script to a file
    with open(output_filename, 'w') as file:
        file.write(slurm_content)

    print(f"SLURM script '{output_filename}' has been created.")

def main():
    parser = argparse.ArgumentParser(description='Generate a SLURM script for fq to bam conversion using Parabricks.')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input fq.gz files')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to output the BAM files')
    parser.add_argument('--output_filename', type=str, default='parabricks_job.slurm', help='Filename for the generated SLURM script')
    parser.add_argument('--read_group_prefix', type=str, default='group', help='Prefix for read group IDs and PUs')

    args = parser.parse_args()

    generate_slurm_script(args.input_dir, args.output_dir, args.output_filename, args.read_group_prefix)

if __name__ == '__main__':
    main()
