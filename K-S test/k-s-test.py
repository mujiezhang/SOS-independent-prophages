
import numpy as np
from scipy import stats

def ks_test_two_samples(list1, list2, alpha=0.05):
    """执行双样本Kolmogorov-Smirnov检验并返回详细结果"""
    arr1 = np.array(list1)
    arr2 = np.array(list2)
    
    ks_stat, p_value = stats.ks_2samp(arr1, arr2)
    
    if p_value < alpha:
        interpretation = f"在{alpha}显著性水平下拒绝原假设: 两组数据分布不同"
    else:
        interpretation = f"在{alpha}显著性水平下无法拒绝原假设: 两组数据分布可能相似"
    
    stats1 = {
        'n': len(arr1),
        'mean': np.mean(arr1),
        'median': np.median(arr1),
        'std': np.std(arr1),
        'min': np.min(arr1),
        'max': np.max(arr1)
    }
    
    stats2 = {
        'n': len(arr2),
        'mean': np.mean(arr2),
        'median': np.median(arr2),
        'std': np.std(arr2),
        'min': np.min(arr2),
        'max': np.max(arr2)
    }
    
    if ks_stat < 0.1:
        effect_size = "可忽略的差异"
    elif ks_stat < 0.2:
        effect_size = "小差异"
    elif ks_stat < 0.3:
        effect_size = "中等差异"
    else:
        effect_size = "大差异"
    
    return {
        'test_name': "Kolmogorov-Smirnov双样本检验",
        'statistic': ks_stat,
        'p_value': p_value,
        'alpha': alpha,
        'interpretation': interpretation,
        'effect_size': effect_size,
        'group1_stats': stats1,
        'group2_stats': stats2
    }

def print_ks_results(results):
    print(f"\n检验名称: {results['test_name']}")
    print(f"K-S统计量 (D值): {results['statistic']:.6f}")
    print(f"P值: {results['p_value']:.6f}")
    print(f"显著性水平(α): {results['alpha']}")
    print(f"效应量: {results['effect_size']}")
    print(f"\n结论: {results['interpretation']}")
    print("=" * 70)

if __name__ == "__main__":
    # 第一个数据集处理
    list1, list2 = [], []
    with open(r"key-gene-data.tsv") as f1:
        for line in f1:
            parts = line.strip().split('\t')
            if parts[1].startswith('SiPs'):
                list1.append(float(parts[-1]))
            elif parts[1].startswith('SdPs'):
                list2.append(float(parts[-1]))
    
    results1 = ks_test_two_samples(list1, list2)
    print_ks_results(results1)

    # 第二个数据集处理
    list3, list4 = [], []
    list5, list6 = [], []
    with open(r"sdp-sip.txt") as f2:
        for line in f2:
            parts = line.strip().split('\t')
            if parts[-1].startswith('SiPs'):
                list3.append(float(parts[2]))
                list5.append(float(parts[1]))
            elif parts[-1].startswith('SdPs'):
                list4.append(float(parts[2]))
                list6.append(float(parts[1]))
    
    results2 = ks_test_two_samples(list3, list4)
    print_ks_results(results2)
    results3 = ks_test_two_samples(list5, list6)
    print_ks_results(results3)
