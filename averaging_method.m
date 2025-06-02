close all;
clear;
clc;

% 数据矩阵
data = [11.6406  10.0445  1.5961  11.7609  10.0445  1.7164  12.8301  10.0445  2.7856
        16.175   15.0431  1.1319  16.3426  15.0431  1.2995  17.7941  15.0431  2.751
        18.8173  18.0426  0.7747  19.0123  18.0426  0.9697  20.7329  18.0426  2.6903
        22.352   22.0589  0.2931  22.5839  22.0589  0.5250  24.5152  22.0589  2.4563
        27.8982  28.0683  0.1701  28.1851  28.0683  0.1168  30.3636  28.0683  2.2953];

% 计算 present method 误差百分比
we = sum(data(:,3)) / sum(data(:,2)) / 5 * 100;
fprintf('Present method: %.4f%%\n', we);

% 计算 Par method 误差百分比
Par = sum(data(:,6)) / sum(data(:,5)) / 5 * 100;
fprintf('Par method: %.4f%%\n', Par);

% 计算 Han method 误差百分比
Han = sum(data(:,9)) / sum(data(:,8)) / 5 * 100;
fprintf('Han method: %.4f%%\n', Han);

% 计算 MAE（Mean Absolute Error）
MAE_present = mean(data(:,3));
MAE_Par = mean(data(:,6));
MAE_Han = mean(data(:,9));

fprintf('MAE (Present method): %.4f\n', MAE_present);
fprintf('MAE (Par method): %.4f\n', MAE_Par);
fprintf('MAE (Han method): %.4f\n', MAE_Han);

% 计算 MASE（Mean Absolute Scaled Error）
% 使用 naive 预测差分作为比例因子
naive_differences = abs(diff(data(:,2))); % 计算相邻真实值的差值
scaling_factor = mean(naive_differences); % 计算缩放因子（均值）

MASE_present = MAE_present / scaling_factor;
MASE_Par = MAE_Par / scaling_factor;
MASE_Han = MAE_Han / scaling_factor;

fprintf('MASE (Present method): %.4f\n', MASE_present);
fprintf('MASE (Par method): %.4f\n', MASE_Par);
fprintf('MASE (Han method): %.4f\n', MASE_Han);



% 相对误差（Relative Error）向量
RelErr_present = data(:,3) ./ abs(data(:,2));
RelErr_Par = data(:,6) ./ abs(data(:,5));
RelErr_Han = data(:,9) ./ abs(data(:,8));

% 平均相对误差（Mean Relative Error, MRE）
MRE_present = mean(RelErr_present);
MRE_Par = mean(RelErr_Par);
MRE_Han = mean(RelErr_Han);

fprintf('Mean Relative Error (Present method): %.4f\n', MRE_present);
fprintf('Mean Relative Error (Par method): %.4f\n', MRE_Par);
fprintf('Mean Relative Error (Han method): %.4f\n', MRE_Han);

