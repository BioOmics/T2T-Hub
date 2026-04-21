import argparse
from Bio import SeqIO

def calculate_gc_content(sequence):
    """
    Calculate the GC content percentage of the sequence, excluding 'N' bases
    :param sequence: Input DNA sequence
    :return: GC content percentage
    """
    # 计算G和C的总数 (calculate the total count of G and C)
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    
    # 计算A、T、C、G的总数，排除'N' (calculate the total count of A, T, C, G, excluding 'N')
    total_bases = sequence.count('A') + sequence.count('T') + g_count + c_count
    
    if total_bases == 0:
        return 0.0  # 防止除零错误 (avoid division by zero)
    
    # 计算GC含量的百分比 (calculate the GC content percentage)
    gc_content = (g_count + c_count) / total_bases * 100
    return gc_content

def calculate_genome_gc(fasta_file):
    """
    Calculate the GC content of the entire genome
    :param fasta_file: Path to the input FASTA file
    :return: GC content percentage of the entire genome
    """
    total_sequence = ""
    
    # 合并所有染色体的序列 (merge the sequences of all chromosomes)
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_sequence += str(record.seq)
    
    # 计算全基因组GC含量 (calculate the GC content of the entire genome)
    return calculate_gc_content(total_sequence)

def calculate_chromosome_gc(fasta_file):
    """
    Calculate the GC content for each chromosome in the FASTA file, excluding 'N' bases
    :param fasta_file: Path to the input FASTA file
    :return: Dictionary with chromosome names as keys and GC content percentages as values
    """
    gc_content_per_chromosome = {}
    
    # 读取FASTA文件中的每条序列 (read each sequence in the FASTA file)
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = record.id
        sequence = str(record.seq)
        
        # 计算该染色体序列的GC含量 (calculate the GC content of the chromosome sequence)
        gc_content = calculate_gc_content(sequence)
        gc_content_per_chromosome[chromosome] = gc_content
    
    return gc_content_per_chromosome

def main():
    # 使用argparse解析命令行参数 (Use argparse to parse command-line arguments)
    parser = argparse.ArgumentParser(description="Calculate GC content for each chromosome in a FASTA file.")
    parser.add_argument("fasta_file", help="Path to the FASTA file containing the genome sequences")
    args = parser.parse_args()
    
    # 获取FASTA文件路径 (Get the path to the FASTA file)
    fasta_file = args.fasta_file
    
    # 计算全基因组GC含量 (Calculate the GC content of the entire genome)
    genome_gc = calculate_genome_gc(fasta_file)
    print(f"Genome\t{genome_gc:.2f}%")
    
    # 计算每个染色体的GC含量 (Calculate the GC content for each chromosome)
    gc_content = calculate_chromosome_gc(fasta_file)
    
    # 输出每个染色体的GC含量 (Output the GC content for each chromosome)
    for chromosome, gc in gc_content.items():
        print(f"{chromosome}\t{gc:.2f}%")

if __name__ == "__main__":
    main()
