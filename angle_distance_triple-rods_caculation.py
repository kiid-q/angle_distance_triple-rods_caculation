import numpy as np
import pandas as pd
from scipy.spatial import KDTree

# ==========================================
# 1. 参数设置区域 (在此修改你的阈值)
# ==========================================
PARAMS = {
    'rod_file': 'rods.txt',              # 杆状蛋白数据
    'lys_file': 'lysosome_membrane.txt',  # 溶酶体膜数据
    'er_file': 'ER_membrane.txt',        # 内质网膜数据
    'output_name': 'analysis_results.csv',# 输出文件名
    
    'r': 15,                             # 膜拟合搜索半径 (影响 Grade)
    'contact_threshold': 2.371,            # 判定“接触”的距离阈值
    'triple_threshold': 13.04             # 判定“三连杆”的中心点间距
}

# ==========================================
# 2. 核心计算逻辑
# ==========================================

def get_fit_grade(linearity):
    """根据线性度判定等级"""
    if linearity > 0.95: return "A"
    elif linearity > 0.85: return "B"
    else: return "C/D"

def analyze_membrane_logic(p_base, p_tip, mem_pts, r):
    """针对单膜进行局部 2D SVD 拟合"""
    # 筛选 Z 平层
    z_layer = int(round(p_base[2]))
    layer_mask = np.abs(mem_pts[:, 2] - z_layer) < 0.5
    layer_pts = mem_pts[layer_mask, :2]
    if len(layer_pts) < 3: return None

    # 局部半径搜索
    tree_2d = KDTree(layer_pts)
    idx = tree_2d.query_ball_point(p_base[:2], r)
    if len(idx) < 3: return None

    # SVD 拟合直线
    local_pts = layer_pts[idx]
    centroid = local_pts.mean(axis=0)
    _, s, vh = np.linalg.svd(local_pts - centroid)
    
    line_vec, norm_vec = vh[0], vh[1]
    v_rod = p_tip[:2] - p_base[:2]
    
    # 计算几何指标
    cos_theta = np.abs(np.dot(v_rod, line_vec)) / (np.linalg.norm(v_rod) + 1e-9)
    angle = np.degrees(np.arccos(np.clip(cos_theta, 0, 1)))
    dist_to_line = np.abs(np.dot(p_base[:2] - centroid, norm_vec))
    linearity = 1.0 - (s[1] / s[0]) if s[0] > 0 else 0
    rmse = s[1] / np.sqrt(len(local_pts))

    return {
        'Angle': round(angle, 2),
        'Dist': round(dist_to_line, 3),
        'RMSE': round(rmse, 3),
        'Linearity': round(linearity, 3),
        'Points': len(local_pts),
        'Grade': get_fit_grade(linearity)
    }

def run_analysis():
    print("正在读取数据并初始化...")
    try:
        rods = pd.read_csv(PARAMS['rod_file'], sep='\s+', names=['X', 'Y', 'Z']).values
        lys_mem = pd.read_csv(PARAMS['lys_file'], sep='\s+', names=['X', 'Y', 'Z']).values
        er_mem = pd.read_csv(PARAMS['er_file'], sep='\s+', names=['X', 'Y', 'Z']).values
    except Exception as e:
        print(f"文件读取错误: {e}")
        return

    results = []
    rod_midpoints = []
    lys_tree_3d = KDTree(lys_mem)
    er_tree_3d = KDTree(er_mem)

    # 遍历每一根杆
    for i in range(0, len(rods) - 1, 2):
        p1, p2 = rods[i], rods[i+1]
        midpt = (p1 + p2) / 2
        rod_midpoints.append(midpt)

        # 确定相对于不同膜的基准点（最近原则）
        d1_l, _ = lys_tree_3d.query(p1); d2_l, _ = lys_tree_3d.query(p2)
        p_base_l, p_tip_l = (p1, p2) if d1_l < d2_l else (p2, p1)
        
        d1_e, _ = er_tree_3d.query(p1); d2_e, _ = er_tree_3d.query(p2)
        p_base_e, p_tip_e = (p1, p2) if d1_e < d2_e else (p2, p1)

        # 执行计算
        lys_res = analyze_membrane_logic(p_base_l, p_tip_l, lys_mem, PARAMS['r'])
        er_res = analyze_membrane_logic(p_base_e, p_tip_e, er_mem, PARAMS['r'])

        row = {'ID': i // 2 + 1, 'Z_Layer': int(round(midpt[2]))}
        
        # 溶酶体数据整合
        if lys_res:
            row.update({f'Lys_{k}': v for k, v in lys_res.items()})
            row['Lys_Contact'] = 1 if lys_res['Dist'] < PARAMS['contact_threshold'] else 0
        else: row['Lys_Contact'] = 0
        
        # 内质网数据整合
        if er_res:
            row.update({f'ER_{k}': v for k, v in er_res.items()})
            row['ER_Contact'] = 1 if er_res['Dist'] < PARAMS['contact_threshold'] else 0
        else: row['ER_Contact'] = 0
        
        results.append(row)

    # 三连杆分析
    rod_tree = KDTree(rod_midpoints)
    for idx, pt in enumerate(rod_midpoints):
        neighbors = rod_tree.query_ball_point(pt, PARAMS['triple_threshold'])
        results[idx]['Is_Triple_Rod'] = 1 if len(neighbors) >= 3 else 0

    # 生成 DataFrame 并计算统计量
    df = pd.DataFrame(results)
    total = len(df)
    lys_ratio = df['Lys_Contact'].sum() / total
    er_ratio = df['ER_Contact'].sum() / total
    triple_ratio = df['Is_Triple_Rod'].sum() / total

    # 构建汇总行
    summary_data = [
        {}, # 空行
        {'ID': '--- STATS SUMMARY ---'},
        {'ID': 'Total Rod Count', 'Z_Layer': total},
        {'ID': 'Lysosome Contact %', 'Z_Layer': f"{lys_ratio:.2%}"},
        {'ID': 'ER Contact %', 'Z_Layer': f"{er_ratio:.2%}"},
        {'ID': 'Triple Rod %', 'Z_Layer': f"{triple_ratio:.2%}"}
    ]
    df_summary = pd.DataFrame(summary_data)
    final_df = pd.concat([df, df_summary], ignore_index=True)

    # 导出
    final_df.to_csv(PARAMS['output_name'], index=False, encoding='utf-8-sig')
    
    print("\n" + "="*30)
    print(f"分析完成！")
    print(f"总计处理杆数: {total}")
    print(f"溶酶体接触比例: {lys_ratio:.2%}")
    print(f"内质网接触比例: {er_ratio:.2%}")
    print(f"三连杆结构比例: {triple_ratio:.2%}")
    print(f"数据已保存至: {PARAMS['output_name']}")
    print("="*30)

if __name__ == "__main__":
    run_analysis()