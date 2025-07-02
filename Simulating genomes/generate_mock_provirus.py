#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模拟不同完整度和污染度的病毒基因组
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


def extract_provirus(host_contigs, provirus_positions):
    """从宿主基因组中提取原病毒序列"""
    provirus_sequences = []
    
    for pos_info in provirus_positions:
        contig_id = pos_info['contig_id']
        start = pos_info['start']
        end = pos_info['end']
        provirus_id = pos_info.get('provirus_id') or f"provirus_{len(provirus_sequences)+1}"
        
        # 查找对应的contig
        target_contig = None
        for contig in host_contigs:
            if contig.id == contig_id or contig_id in contig.id:
                target_contig = contig
                break
        
        if target_contig is None:
            print(f"警告: 找不到contig {contig_id}")
            continue
        
        # 提取原病毒序列
        contig_seq = str(target_contig.seq)
        if end > len(contig_seq):
            print(f"警告: 原病毒位置超出contig长度 ({contig_id}: {end} > {len(contig_seq)})")
            end = len(contig_seq)
        
        if start >= end:
            print(f"警告: 原病毒位置无效 ({contig_id}: {start} >= {end})")
            continue
        
        provirus_seq = contig_seq[start-1:end]  # 转换为0-based索引
        provirus_record = SeqRecord(
            Seq(provirus_seq),
            id=provirus_id,
            description=f"provirus_from_{contig_id}_{start}_{end}"
        )
        provirus_sequences.append(provirus_record)
        
        print(f"提取原病毒: {provirus_id} (长度: {len(provirus_seq):,} bp)")
    
    return provirus_sequences


def simulate_provirus_completeness(provirus_sequences, completeness):
    """通过保留连续子片段模拟原病毒完整度"""
    incomplete_sequences = []
    
    for provirus in provirus_sequences:
        original_seq = str(provirus.seq)
        original_len = len(original_seq)
        target_len = int(original_len * completeness)
        
        if target_len >= original_len:
            # 完整度100%或更高，保持原序列
            incomplete_seq = original_seq
            kept_start = 1
            kept_end = original_len
        else:
            # 随机选择连续子片段的起始位置
            max_start = original_len - target_len
            start_pos = random.randint(0, max_start)
            end_pos = start_pos + target_len
            
            incomplete_seq = original_seq[start_pos:end_pos]
            kept_start = start_pos + 1  # 转换为1-based坐标用于描述
            kept_end = end_pos
        
        incomplete_record = SeqRecord(
            Seq(incomplete_seq),
            id=f"{provirus.id}_incomplete",
            description=f"{provirus.description}_kept_{kept_start}_{kept_end}_completeness_{completeness:.2f}"
        )
        incomplete_sequences.append(incomplete_record)
        
        actual_completeness = len(incomplete_seq) / original_len
        print(f"原病毒 {provirus.id}: 原长度 {original_len:,} bp -> 保留片段 {len(incomplete_seq):,} bp (位置: {kept_start}-{kept_end}, 完整度: {actual_completeness:.1%})")
    
    return incomplete_sequences


def add_host_flanking_regions(provirus_sequences, host_contigs, contamination_rate, fragment_size_range=(1000, 10000)):
    """在原病毒两端添加宿主基因组片段，组成连续序列"""
    if contamination_rate == 0:
        return provirus_sequences
    
    final_sequences = []
    
    for provirus in provirus_sequences:
        provirus_seq = str(provirus.seq)
        provirus_length = len(provirus_seq)
        
        # 计算宿主片段总长度
        total_host_length = int(provirus_length * contamination_rate / (1 - contamination_rate))
        
        # 决定宿主片段的分布模式
        flank_mode = random.choice(['both', 'left', 'right'])
        
        left_host_seq = ""
        right_host_seq = ""
        
        if flank_mode == 'both':
            # 两端都有宿主片段
            left_length = total_host_length // 2
            right_length = total_host_length - left_length
            
            if left_length > 0:
                left_host_seq = get_random_host_fragment(host_contigs, left_length, fragment_size_range)
            if right_length > 0:
                right_host_seq = get_random_host_fragment(host_contigs, right_length, fragment_size_range)
                
        elif flank_mode == 'left':
            # 只有左端有宿主片段
            left_host_seq = get_random_host_fragment(host_contigs, total_host_length, fragment_size_range)
            
        else:  # right
            # 只有右端有宿主片段
            right_host_seq = get_random_host_fragment(host_contigs, total_host_length, fragment_size_range)
        
        # 组装最终序列：[左侧宿主] + [原病毒] + [右侧宿主]
        final_seq = left_host_seq + provirus_seq + right_host_seq
        
        # 计算实际污染度
        actual_host_length = len(left_host_seq) + len(right_host_seq)
        actual_cont_rate = actual_host_length / len(final_seq) if len(final_seq) > 0 else 0
        
        # 创建最终序列记录
        final_record = SeqRecord(
            Seq(final_seq),
            id=f"{provirus.id}_with_flanks",
            description=f"{provirus.description}_flanked_left{len(left_host_seq)}_right{len(right_host_seq)}_cont{actual_cont_rate:.3f}"
        )
        
        final_sequences.append(final_record)
        
        print(f"原病毒 {provirus.id}:")
        print(f"  - 原病毒长度: {provirus_length:,} bp")
        print(f"  - 左侧宿主: {len(left_host_seq):,} bp")
        print(f"  - 右侧宿主: {len(right_host_seq):,} bp")
        print(f"  - 最终长度: {len(final_seq):,} bp")
        print(f"  - 实际污染度: {actual_cont_rate:.1%}")
    
    return final_sequences


def get_random_host_fragment(host_contigs, target_length, fragment_size_range):
    """从宿主基因组中获取指定长度的随机片段"""
    if target_length <= 0:
        return ""
    
    host_seq = ""
    remaining_length = target_length
    attempts = 0
    
    while remaining_length > 0 and attempts < 100:  # 防止无限循环
        # 随机选择一个contig
        contig = random.choice(host_contigs)
        contig_seq = str(contig.seq)
        contig_len = len(contig_seq)
        
        if contig_len < 100:  # 跳过太短的contig
            attempts += 1
            continue
        
        # 确定这次提取的片段长度
        max_extract = min(fragment_size_range[1], contig_len, remaining_length)
        min_extract = min(fragment_size_range[0], max_extract)
        
        if min_extract <= 0:
            break
            
        fragment_length = random.randint(min_extract, max_extract)
        
        # 随机选择起始位置
        start_pos = random.randint(0, contig_len - fragment_length)
        fragment = contig_seq[start_pos:start_pos + fragment_length]
        
        host_seq += fragment
        remaining_length -= fragment_length
        attempts += 1
    
    return host_seq


def save_sequences(sequences, output_path):
    """保存序列到文件"""
    try:
        SeqIO.write(sequences, output_path, "fasta")
        total_length = sum(len(seq.seq) for seq in sequences)
        print(f"保存到: {output_path}")
        print(f"包含 {len(sequences)} 个序列, 总长度 {total_length:,} bp")
        return True
    except Exception as e:
        print(f"保存失败: {e}")
        return False


def parse_positions(position_str):
    """解析原病毒位置信息"""
    positions = []
    
    # 支持多种格式：
    # 1. "contig1:1000-2000"
    # 2. "contig1:1000-2000,contig2:3000-4000"
    # 3. "contig1:1000-2000:provirus1,contig2:3000-4000:provirus2"
    
    for pos_part in position_str.split(','):
        pos_part = pos_part.strip()
        
        if ':' in pos_part:
            parts = pos_part.split(':')
            contig_id = parts[0]
            pos_range = parts[1]
            provirus_id = parts[2] if len(parts) > 2 else None
            
            if '-' in pos_range:
                start, end = map(int, pos_range.split('-'))
                positions.append({
                    'contig_id': contig_id,
                    'start': start,
                    'end': end,
                    'provirus_id': provirus_id
                })
        else:
            print(f"警告: 无法解析位置信息: {pos_part}")
    
    return positions


def simulate_provirus_genome(host_genome, provirus_positions, completeness, contamination=0, 
                           output=None, seed=None):
    """主要的原病毒模拟函数"""
    if seed is not None:
        random.seed(seed)
        print(f"随机种子: {seed}")
    
    print("=" * 60)
    print("开始原病毒基因组模拟")
    print("=" * 60)
    
    # 1. 读取宿主基因组
    print("\n1. 读取宿主基因组...")
    host_contigs = read_genome(host_genome)
    if not host_contigs:
        return False
    
    # 2. 解析原病毒位置
    print("\n2. 解析原病毒位置...")
    if isinstance(provirus_positions, str):
        positions = parse_positions(provirus_positions)
    else:
        positions = provirus_positions
    
    print(f"找到 {len(positions)} 个原病毒位置")
    for pos in positions:
        print(f"  - {pos['contig_id']}:{pos['start']}-{pos['end']}")
    
    # 3. 提取原病毒序列
    print("\n3. 提取原病毒序列...")
    provirus_seqs = extract_provirus(host_contigs, positions)
    if not provirus_seqs:
        print("错误: 没有成功提取原病毒序列")
        return False
    
    # 4. 模拟完整度（保留连续子片段）
    print(f"\n4. 模拟完整度 ({completeness:.1%}) - 保留连续子片段...")
    incomplete_seqs = simulate_provirus_completeness(provirus_seqs, completeness)
    
    # 5. 添加宿主侧翼区域
    if contamination > 0:
        print(f"\n5. 添加宿主侧翼区域 ({contamination:.1%})...")
        final_seqs = add_host_flanking_regions(incomplete_seqs, host_contigs, contamination)
    else:
        final_seqs = incomplete_seqs
    
    # 6. 保存结果
    if output is None:
        base_name = os.path.splitext(os.path.basename(host_genome))[0]
        output = f"{base_name}_provirus_c{completeness:.2f}_cont{contamination:.2f}.fasta"
    
    print(f"\n6. 保存结果...")
    success = save_sequences(final_seqs, output)
    
    if success:
        print("\n" + "=" * 60)
        print("原病毒模拟完成!")
        print("=" * 60)
    
    return success


def main():
    parser = argparse.ArgumentParser(description="script for simulating proviral genomes")
    
    parser.add_argument("host_genome", help="宿主基因组文件 (FASTA)")
    parser.add_argument("-p", "--positions", required=True,
                       help="原病毒位置 (格式: contig:start-end 或 contig:start-end:id)")
    parser.add_argument("-c", "--completeness", type=float, required=True,
                       help="完整度 (0-1)")
    parser.add_argument("-cont", "--contamination", type=float, default=0,
                       help="污染度 (0-1)")
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
    
    if not os.path.exists(args.host_genome):
        print(f"错误: 宿主基因组文件不存在: {args.host_genome}")
        sys.exit(1)
    
    # 运行模拟
    success = simulate_provirus_genome(
        host_genome=args.host_genome,
        provirus_positions=args.positions,
        completeness=args.completeness,
        contamination=args.contamination,
        output=args.output,
        seed=args.seed
    )
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
