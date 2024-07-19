import os
import argparse

def generate_slurm_script(input_dir, output_dir, output_filename):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Find all .bam files in the input directory
    bam_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.bam')])

    # Generate commands for each bam file
    commands = []
    for bam_file in bam_files:
        variant_path = os.path.join(output_dir, f"{os.path.splitext(bam_file)[0]}.vcf")
        command = f"""
singularity run --cleanenv --nv \\
    --bind /home/jl2791:/home/jl2791,/scratch/jl2791:/scratch/jl2791,/projectsp/f_ak1833_1:/projectsp/f_ak1833_1,/projects:/projects,/cache:/cache \\
    /projects/community/singularity.images/nvidia/clara-parabricks_4.3.1-1.sif pbrun haplotypecaller \\
    --ref /projectsp/f_ak1833_1/Reference_Data/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \\
    --in-bam {os.path.join(input_dir, bam_file)} \\
    --out-variants {variant_path} \\
    --num-gpus 4 --run-partition --htvc-low-memory --num-htvc-threads 10
"""
        commands.append(command)

    # Create the SLURM script content
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=parabricks_haplotypecaller
#SBATCH --partition=gpu
#SBATCH --output={output_dir}/parabricks_haplotypecaller_%j.out
#SBATCH --error={output_dir}/parabricks_haplotypecaller_%j.err
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
{''.join(commands)}
    """

    # Write the script to a file
    with open(output_filename, 'w') as file:
        file.write(slurm_content)

    print(f"SLURM script '{output_filename}' has been created.")

def main():
    parser = argparse.ArgumentParser(description='Generate a SLURM script for running HaplotypeCaller using Parabricks.')
    parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input BAM files')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to output the VCF files')
    parser.add_argument('--output_filename', type=str, default='02.parabricks_haplotypecaller_job.slurm', help='Filename for the generated SLURM script')

    args = parser.parse_args()

    generate_slurm_script(args.input_dir, args.output_dir, args.output_filename)

if __name__ == '__main__':
    main()
