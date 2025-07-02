#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模拟不同完整度和污染度的细菌基因组
"""

import os
import sys
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_genome(genome_path):
    """读取基因组文件"""
    try:
        records = list(SeqIO.parse(genome_path, "fasta"))
        print(f"读取基因组: {genome_path}, 包含 {len(records)} 个contig")
        total_length = sum(len(record.seq) for record in records)
        print(f"基因组总长度: {total_length:,} bp")
        return records
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None


def fragment_genome(contigs, min_size=5000, max_size=50000):
    """将基因组打碎成片段"""
    fragments = []
    fragment_id = 1
    
    for contig in contigs:
        contig_seq = str(contig.seq)
        contig_len = len(contig_seq)
        
        # 短contig直接保留
        if contig_len <= max_size:
            fragment = SeqRecord(
                Seq(contig_seq),
                id=f"fragment_{fragment_id:06d}",
                description=f"from_{contig.id}"
            )
            fragments.append(fragment)
            fragment_id += 1
            continue
        
        # 长contig打碎
        start = 0
        while start < contig_len:
            fragment_size = random.randint(min_size, max_size)
            end = min(start + fragment_size, contig_len)
            
            # 确保最后片段不会太小
            if contig_len - end < min_size and end < contig_len:
                end = contig_len
            
            if end - start >= min_size:
                fragment_seq = contig_seq[start:end]
                fragment = SeqRecord(
                    Seq(fragment_seq),
                    id=f"fragment_{fragment_id:06d}",
                    description=f"from_{contig.id}_{start}_{end}"
                )
                fragments.append(fragment)
                fragment_id += 1
            
            start = end
    
    print(f"基因组打碎成 {len(fragments)} 个片段")
    return fragments


def select_by_completeness(fragments, completeness):
    """根据完整度选择片段"""
    total_length = sum(len(f.seq) for f in fragments)
    target_length = int(total_length * completeness)
    
    # 随机打乱并选择
    shuffled = fragments.copy()
    random.shuffle(shuffled)
    
    selected = []
    current_length = 0
    
    for fragment in shuffled:
        if current_length >= target_length:
            break
        selected.append(fragment)
        current_length += len(fragment.seq)
    
    actual_completeness = current_length / total_length
    print(f"完整度: 目标 {completeness:.1%}, 实际 {actual_completeness:.1%}")
    print(f"选择了 {len(selected)} 个片段, 总长度 {current_length:,} bp")
    
    return selected


def add_contamination(main_fragments, contaminant_files, contamination_rate):
    """添加污染片段"""
    if contamination_rate == 0 or not contaminant_files:
        return main_fragments
    
    main_length = sum(len(f.seq) for f in main_fragments)
    target_cont_length = int(main_length * contamination_rate / (1 - contamination_rate))
    
    # 收集所有污染片段
    all_cont_fragments = []
    for cont_file in contaminant_files:
        if os.path.exists(cont_file):
            cont_contigs = read_genome(cont_file)
            if cont_contigs:
                cont_fragments = fragment_genome(cont_contigs)
                all_cont_fragments.extend(cont_fragments)
    
    if not all_cont_fragments:
        print("警告: 没有找到污染基因组片段")
        return main_fragments
    
    # 随机选择污染片段
    random.shuffle(all_cont_fragments)
    selected_cont = []
    current_cont_length = 0
    
    for fragment in all_cont_fragments:
        if current_cont_length >= target_cont_length:
            break
        fragment.id = f"contaminant_{len(selected_cont)+1:06d}"
        selected_cont.append(fragment)
        current_cont_length += len(fragment.seq)
    
    # 计算实际污染度
    total_length = main_length + current_cont_length
    actual_cont_rate = current_cont_length / total_length
    print(f"污染度: 目标 {contamination_rate:.1%}, 实际 {actual_cont_rate:.1%}")
    print(f"添加了 {len(selected_cont)} 个污染片段")
    
    # 合并并打乱
    all_fragments = main_fragments + selected_cont
    random.shuffle(all_fragments)
    
    return all_fragments


def save_genome(fragments, output_path):
    """保存基因组到文件"""
    try:
        SeqIO.write(fragments, output_path, "fasta")
        total_length = sum(len(f.seq) for f in fragments)
        print(f"保存到: {output_path}")
        print(f"包含 {len(fragments)} 个contig, 总长度 {total_length:,} bp")
        return True
    except Exception as e:
        print(f"保存失败: {e}")
        return False


def simulate_genome(main_genome, completeness, contamination=0, contaminant_genomes=None, 
                   output=None, seed=None):
    """主要模拟函数"""
    if seed is not None:
        random.seed(seed)
        print(f"随机种子: {seed}")
    
    print("=" * 50)
    print("开始基因组模拟")
    print("=" * 50)
    
    # 1. 读取主基因组
    print("\n1. 读取主基因组...")
    main_contigs = read_genome(main_genome)
    if not main_contigs:
        return False
    
    # 2. 打碎基因组
    print("\n2. 打碎基因组...")
    fragments = fragment_genome(main_contigs)
    
    # 3. 根据完整度选择
    print(f"\n3. 选择片段 (完整度: {completeness:.1%})...")
    selected_fragments = select_by_completeness(fragments, completeness)
    
    # 4. 添加污染
    if contamination > 0:
        print(f"\n4. 添加污染 (污染度: {contamination:.1%})...")
        final_fragments = add_contamination(selected_fragments, contaminant_genomes or [], contamination)
    else:
        final_fragments = selected_fragments
    
    # 5. 保存结果
    if output is None:
        base_name = os.path.splitext(os.path.basename(main_genome))[0]
        output = f"{base_name}_mock_c{completeness:.2f}_cont{contamination:.2f}.fasta"
    
    print(f"\n5. 保存结果...")
    success = save_genome(final_fragments, output)
    
    if success:
        print("\n" + "=" * 50)
        print("模拟完成!")
        print("=" * 50)
    
    return success


def main():
    parser = argparse.ArgumentParser(description="scirpt for simulating bacterial genomes")
    
    parser.add_argument("genome", help="主基因组文件 (FASTA)")
    parser.add_argument("-c", "--completeness", type=float, required=True, 
                       help="完整度 (0-1)")
    parser.add_argument("-cont", "--contamination", type=float, default=0,
                       help="污染度 (0-1)")
    parser.add_argument("-cg", "--contaminant-genomes", nargs="+", default=[],
                       help="污染基因组文件列表")
    parser.add_argument("-o", "--output", help="输出文件")
    parser.add_argument("-s", "--seed", type=int, help="随机种子")
    
    args = parser.parse_args()
    
    # 检查参数
    if not 0 < args.completeness <= 1:
        print("错误: 完整度必须在0-1之间")
        sys.exit(1)
    
    if not 0 <= args.contamination < 1:
        print("错误: 污染度必须在0-1之间")
        sys.exit(1)
    
    if not os.path.exists(args.genome):
        print(f"错误: 基因组文件不存在: {args.genome}")
        sys.exit(1)
    
    if args.contamination > 0 and not args.contaminant_genomes:
        print("错误: 指定了污染度但没有污染基因组文件")
        sys.exit(1)
    
    # 运行模拟
    success = simulate_genome(
        main_genome=args.genome,
        completeness=args.completeness,
        contamination=args.contamination,
        contaminant_genomes=args.contaminant_genomes,
        output=args.output,
        seed=args.seed
    )
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
