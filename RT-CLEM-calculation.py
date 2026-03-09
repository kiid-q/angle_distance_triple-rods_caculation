import os
import re
import csv

def export_sums_to_csv(folder_path):
    # 结果存储
    data_list = []
    total_b = 0  # 替换为 B
    total_d = 0  # 替换为 D
    total_i = 0  # 替换为 I
    
    # 【核心修改点】正则表达式：匹配文件名中的 B, D, I 后面的数字
    # 模式变为：_B数字_D数字_I数字
    pattern = re.compile(r'_B(\d+)_D(\d+)_I(\d+)')

    # 遍历文件夹
    filenames = [f for f in os.listdir(folder_path) if f.lower().endswith(".tif")]
    
    for filename in filenames:
        match = pattern.search(filename)
        if match:
            b_val = int(match.group(1))
            d_val = int(match.group(2))
            i_val = int(match.group(3))
            
            # 累加
            total_b += b_val
            total_d += d_val
            total_i += i_val
            
            # 记录单行数据
            data_list.append([filename, b_val, d_val, i_val])

    if not data_list:
        print("错误：未发现匹配 '_B#_D#_I#' 格式的 .tif 文件。")
        return

    # 写入 CSV 文件
    output_file = "image_data_summary_BDI.csv"
    
    try:
        with open(output_file, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            # 写入表头
            writer.writerow(['FileName', 'B_Value', 'D_Value', 'I_Value'])
            # 写入数据
            writer.writerows(data_list)
            # 写入总计行
            writer.writerow(['TOTAL_SUM', total_b, total_d, total_i])
            
        print("-" * 30)
        print(f"成功！CSV 文件已生成: {output_file}")
        print(f"总计结果 -> B:{total_b}, D:{total_d}, I:{total_i}")
        print("-" * 30)
    except PermissionError:
        print("错误：无法写入文件。请先关闭 Excel 中的 image_data_summary_BDI.csv 再运行。")

if __name__ == "__main__":
    export_sums_to_csv(os.getcwd())