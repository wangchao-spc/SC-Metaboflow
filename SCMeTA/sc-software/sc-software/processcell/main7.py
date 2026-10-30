import matplotlib.pyplot as plt
from fisher_py import RawFile
import numpy as np
from fisher_py.data.business import TraceType
from scipy.interpolate import interp1d
from processData import find_peaks_and_valleys, clean_data, find_valley_before_peak_rt, exclude_isotopes, segmented_normalization
import os

# 检查文件夹是否存在，如果不存在则创建
output_folder = "removebgr-0408"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 原始文件夹路径
raw_folder_path = "./data/neg/PFOS-0408"

# 获取原始文件夹下所有的.RAW文件路径
raw_file_paths = [os.path.join(raw_folder_path, filename) for filename in os.listdir(raw_folder_path) if filename.endswith(".RAW")]

target_mass = 885.5524
mass_tolerance_ppm = 5
for raw_file_path in raw_file_paths:
    # 读文件
    raw_file = RawFile(raw_file_path)

    # 提取不包含扩展名的原始文件名
    raw_filename_without_extension = os.path.splitext(os.path.basename(raw_file_path))[0]

    # 读原始数据
    rt, i = raw_file.get_chromatogram(target_mass, mass_tolerance_ppm, TraceType.MassRange)

    plt.figure()
    plt.plot(rt, i)


    # 对原始数据做清洗
    cleand_i = clean_data(i)
    plt.figure()
    plt.plot(rt, cleand_i)
    print('max_cleand_i:',np.max(cleand_i))
    

   # 使用新的分段归一化方法替换原有的归一化代码
    normalized_ri = segmented_normalization(rt, cleand_i, segment_length=5)
    plt.figure()
    plt.plot(rt, normalized_ri)
    print(len(rt), rt)
    print('max_normalized_ri:', np.max(normalized_ri))  # 替换为normalized_ri的最大值


    # 找到清洗后的归一化的数据的所有的峰值和谷值
    peak_indexes,valley_indexes = find_peaks_and_valleys(rt, normalized_ri, threshold=0.1, order=1)
    print('Peaks:', peak_indexes)
    print('Valleys:', valley_indexes)

    # 找到所有的峰值，和它对应的前置谷值
    valleyPeakRtTuple = find_valley_before_peak_rt(peak_indexes,valley_indexes,rt)
    print('valleyPeakTuple:',valleyPeakRtTuple)

    # 在遍历峰值之前，初始化一个列表来收集所有峰的质心数据
    all_centroids = []

    # 遍历所有的峰值
    for i, valleyPeakRtItem in enumerate(valleyPeakRtTuple):
    
        # 获取当前峰值的质谱图数据
        mz_peak, i2_peak, charges_peak, real_rt_peak = raw_file.get_scan_ms1(valleyPeakRtItem[1])
    
        # 获取当前谷值的质谱图数据
        mz_valley, i2_valley, charges_valley, real_rt_valley = raw_file.get_scan_ms1(valleyPeakRtItem[3])
    
        print('len mz_peak mz_valley',len(mz_peak),len(mz_valley))
        plt.figure()
    
        # 因为峰值的质谱图和谷值的质谱图的点数不一样，所以需要做插值处理
        x1 = np.array([])
        y1 = np.array([])
        x2 = np.array([])
        y2 = np.array([])
        y_diff = np.array([])
        y_diff[y_diff < 0] = 0  # 将差值小于0的数值归为0
   
    
        # 插值处理，点数少的做插值
        if len(mz_valley) < len(mz_peak):
            print('len mz_valley < len mz_peak')
            x1 = mz_peak
            y1 = i2_peak
            x2 = mz_valley # 点少，需要插值
            y2 = i2_valley # 点少，需要插值
        
            # 插值用的是slinear，一阶样条插值
            f = interp1d(x2, y2, kind='slinear')
            x1_clipped = np.clip(x1, np.min(x2), np.max(x2))
            y2_interp = f(x1_clipped)
            plt.plot(x1_clipped, y2_interp, label='y2_interp')
        
            # 计算两个折线图的差值
            y_diff =  y1 - y2_interp
        
            # 将差值小于0的数值替换为0
            y_diff = np.where(y_diff < 0, 0, y_diff)


        else:
            print('len mz_valley > len mz_peak')
            x1 = mz_valley
            y1 = i2_valley
            x2 = mz_peak # 点少，需要插值
            y2 = i2_peak # 点少，需要插值
        
            # 插值用的是slinear，一阶样条插值
            f = interp1d(x2, y2, kind='slinear')                                                                    
            x1_clipped = np.clip(x1, np.min(x2), np.max(x2))
            y2_interp = f(x1_clipped)
            plt.plot(x1_clipped, y2_interp, label='y2_interp')
        
            # 计算两个折线图的差值
            y_diff = y2_interp - y1
      
            y_diff = np.where(y_diff < 0, 0, y_diff)


        plt.plot(x1, y1, label='y1')
        plt.plot(x2, y2, label='y2')    
        plt.plot(x1, y_diff, label='diff')
        plt.legend()
        # 这部分代码处理了质谱图数据的插值和差值计算，接下来是将结果保存到CSV文件
        y_diff_max = np.max(y_diff)  # 计算y_diff的最大值
        y_diff_percent = (y_diff / y_diff_max) * 100  # 计算相对百分含量
    
        # 将处理后的数据（质量/强度对）存储到一个数组中
        peaks_data = np.column_stack((x1, y_diff))
        peaks_data = peaks_data[peaks_data[:, 1] > 0]  # 去除强度小于等于0的数据

        # 直接使用peaks_data中的数据作为质心数据，而不进行同位素排除
        centroids = [(mz, intensity) for mz, intensity in peaks_data]
        all_centroids.extend(centroids)
        # 计算百分比
        y_diff_max = np.max(peaks_data[:, 1])  # 找到最大强度值
        percentages = (peaks_data[:, 1] / y_diff_max) * 100  # 计算百分比
    
            
        # 每个峰值生成一个CSV文件
        csv_filename = os.path.join(output_folder, f"{raw_filename_without_extension}_peak_{i+1}.csv")
        with open(csv_filename, 'w') as f:
            f.write("m/z,Intensity,Percentage\n")
            for (mz, intensity), percentage in zip(centroids, percentages):
                f.write(f"{mz},{intensity},{percentage:.5f}\n")

 
