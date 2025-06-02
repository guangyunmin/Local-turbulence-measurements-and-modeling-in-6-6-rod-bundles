clc; clear; close all;

% **读取并处理流场数据**
% 你的数据点
data = [
0	0	0.87648
0	8.3333	0.87226
0	25	0.81886
0	16.66667	0.83717
-8.33333	16.66667	0.8246
-16.66667	16.66667	0.8227
-16.6667	25	0.80727
8.33333	16.66667	0.79523
16.66667	16.66667	0.80331
16.6667	25	0.7982
16.6667	41.6667	0.76611
0	33.33333	0.8151
-8.33333	33.33333	0.82813
-16.66667	33.33333	0.8242
-25	33.33333	0.80206
-33.33333	33.33333	0.81298
8.33333	33.33333	0.82522
16.66667	33.33333	0.8117
25	33.33333	0.79087
33.33333	33.33333	0.81046
33.3333	41.6667	0.76817
0	41.6667	0.80085
-16.6667	41.6667	0.76694
-33.3333	41.6667	0.78827
0	48.335	0.81935
-8.33333	48.335	0.74486
-16.66667	48.335	0.80242
-25	48.335	0.77712
-33.33333	48.335	0.80957
-41.66667	48.335	0.73651
8.33333	48.335	0.76335
16.66667	48.335	0.79438
25	48.335	0.78357
33.33333	48.335	0.79003
41.66667	48.335	0.73133
-48.335	48.335	0.74565
48.335	48.335	0.75844
];
save('data3.mat','data')
% =========================
% 🎯 Calculate Subchannel Area
% =========================
P = 16.7;
D = 10;
W = 1.665; 
Across = (100^2-(36*pi*(D/2)^2))/4;

% A-1 to A-25 (1/8 center subchannel, 1 <= i <= 25)
Ac_1_25 = 0.5 * sqrt(5) * P / 4 * sqrt(2) * P / 2 * ...
          sin(pi / 4 - atan(1/2)) - ...
          0.5 * (sqrt(5) * P / 4)^2 * (pi / 4  - atan(1/2));
Asubc_1_25 = 0.5 * (P / 2)^2 - 0.5 * (D / 2)^2 * pi / 4;
Cc_1_25 = Ac_1_25/Asubc_1_25;
Ce_1_25 =1 - Cc_1_25;

% A-26 to A-30 (Type A 1/4 wall subchannel, 26 <= i <= 30)
Ac_26_30 = 0.5 * sqrt((P / 2)^2 + (P / 2 - W)^2) * ...
           sqrt((P / 2)^2 + (P/4 - W / 2)^2) * ...
           sin(atan((P / 2 - W) / (P/2)) - atan((P / 4 - W / 2) / (P / 2))) - ...
           0.5 * ((P / 2)^2 + (P/4 - W/2)^2) * ...
           (atan((P / 2 - W) / (P / 2)) - ...
           atan((P / 4 - W / 2) / (P / 2)));
Asubc_26_30 = 0.5 * (P / 2) * (P / 2 - W) - ...
              0.5 * (D / 2)^2 * atan((P / 2 - W) / (P / 2 ));
Cc_26_30 = Ac_26_30 / Asubc_26_30;
Ce_26_30 =1 - Cc_26_30;

% A-31 to A-35 (Type B 1/4 wall subchannel, 31 <= i <= 35)
Ac_31_35 = 0.5 * sqrt((P/2)^2 + (P/2 - W)^2) * ...
           sqrt((P/4)^2 + (P/2 - W)^2) * ...
           sin(atan(P/2 / (P/2 - W)) - atan(P/4 / (P/2 - W))) - ...
           0.5 * ((P/4)^2 + (P/2 - W)^2) * (atan(P/2 / (P/2 - W)) - ...
           atan(P/4 / (P/2 - W))) + P * W / 8;
Asubc_31_35 = 0.5 * (P/2) * (P/2 + W) - ...
              0.5 * (D / 2)^2 * atan(P/2 / (P/2 - W));
Cc_31_35 = Ac_31_35/Asubc_31_35;
Ce_31_35 =1 - Cc_31_35;
          
% A-36 (1/2 corner subchannel, i = 36)
Ac_36 = 0.5 * (sqrt(5) / 2 * (P / 2 - W)) * ...
        (sqrt(2) * (P / 2 - W)) * ...
        sin(pi/4 - atan(1/2)) - ...
        0.5 * (sqrt(5)/ 2 * (P / 2 - W))^2 * (pi/4 - atan(1/2)) + ...
        0.5 * (P / 2 - W / 2) * W / 2;
Asubc_36 = 0.5 * (P / 2)^2 - 0.5 * (D / 2)^2 * pi / 4;

Cc_36 = Ac_36 / Asubc_36;
Ce_36 =1 - Cc_36;

% Display results
fprintf('Ac_1_25: %.4f\n', Cc_1_25);
fprintf('Asubc_1_25: %.4f\n', Ce_1_25);

fprintf('Ac_26_30: %.4f\n', Cc_26_30);
fprintf('Asubc_26_30: %.4f\n', Ce_26_30);

fprintf('Ac_31_35: %.4f\n', Cc_31_35);
fprintf('Asubc_31_35: %.4f\n',Ce_31_35);

fprintf('Ac_36: %.4f\n', Cc_36);
fprintf('Asubc_36: %.4f\n', Ce_36);


% 提取坐标和值
x = data(:,1);
y = data(:,2);
v = data(:,3);

% 生成所有对称点
x_full = [x; -x; x; -x; y; -y; y; -y];
y_full = [y; y; -y; -y; x; x; -x; -x];
v_full = [v; v; v; v; v; v; v; v];

% **去重**
[unique_points, ~, idx] = unique([x_full, y_full], 'rows');
v_avg = accumarray(idx, v_full, [], @mean);

% 规则网格
[Xq, Yq] = meshgrid(linspace(-50, 50, 100), linspace(-50, 50, 100));

% **插值流场数据**
Vq = griddata(unique_points(:,1), unique_points(:,2), v_avg, Xq, Yq, 'cubic');

% **获取流场的最小 Z 值**
Z_min = min(Vq(:));  % 圆柱底部应从流场的最小 Z 开始
H = 0.3;             % 圆柱高度 0.3


figure;
scatter(x, y, 50, v, 'filled');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('排序前的1/4流场');

    for j = 1:length(x)
        text(x(j), y(j), num2str(v(j), '%.2f'), 'FontSize', 8, 'Color', 'k', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
grid on;
view(-90, 90);  % 旋转 90 度，使 XY 平面对齐
% axis equal;
colormap('jet')


% **步骤 1: 按 Y 从大到小，再按X从大到小排序**
sortedData = sortrows([data(:,1), data(:,2), data(:,3)], [-2, -1]);

% **提取排序后的 X, Y, V**
x_sorted = sortedData(:,1);
y_sorted = sortedData(:,2);
v_sorted = sortedData(:,3);

% **绘制排序后的散点图**
figure;
scatter( x_sorted,y_sorted, 50, v_sorted, 'filled');
colorbar;
xlabel('Y (mm)');
ylabel('X (mm)');
title('排序后的1/4流场');
grid on;
view(-90, 90);  % 旋转 90 度，使 XY 平面对齐
colormap('jet');

% **在图上标注数值**
for i = 1:length(x_sorted)
    text(x_sorted(i), y_sorted(i),sprintf('%.2f', v_sorted(i)), 'FontSize', 8, 'Color', 'k', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

V_corner =  0.657 * v_sorted(1) +  0.223 * v_sorted(2) + ...
    0.657 * v_sorted(13) +  0.223 * v_sorted(12);
V_edge_A = 0.636 * v_sorted(3) +  0.234 * v_sorted(2) + ...
    0.636 * v_sorted(3) +  0.234 * v_sorted(4) + ...
    0.636 * v_sorted(5) +  0.234 * v_sorted(4) + ...
    0.636 * v_sorted(5) +  0.234 * v_sorted(6) + ...
    0.636 * v_sorted(7) +  0.234 * v_sorted(6) + ...
    0.636 * v_sorted(7) +  0.234 * v_sorted(8) + ...
    0.636 * v_sorted(9) +  0.234 * v_sorted(8) + ...
    0.636 * v_sorted(9) +  0.234 * v_sorted(10) + ...
    0.636 * v_sorted(11) +  0.234 * v_sorted(10) + ...
    0.636 * v_sorted(11) +  0.234 * v_sorted(12);



V_edge_B = 2 * (0.659 * v_sorted(3) +  0.222 * v_sorted(14) + ...
    0.659 * v_sorted(5) +  0.222 * v_sorted(15) + ...
    0.659 * v_sorted(7) +  0.222 * v_sorted(16) + ...
    0.659 * v_sorted(9) +  0.222 * v_sorted(17) + ...
    0.659 * v_sorted(11) +  0.222 * v_sorted(18));


V_avg_three_colom = 2*( 0.665 * v_sorted(19) +  0.243 * v_sorted(14) + ...
0.665 * v_sorted(19) +  0.243 * v_sorted(20) + ...
0.665 * v_sorted(21) +  0.243 * v_sorted(20) + ...
0.665 * v_sorted(21) +  0.243 * v_sorted(15) + ...
0.665 * v_sorted(21) +  0.243 * v_sorted(22) + ...
0.665 * v_sorted(23) +  0.243 * v_sorted(22) + ...
0.665 * v_sorted(23) +  0.243 * v_sorted(16) + ...
0.665 * v_sorted(23) +  0.243 * v_sorted(24) + ...
0.665 * v_sorted(25) +  0.243 * v_sorted(24) + ...
0.665 * v_sorted(25) +  0.243 * v_sorted(17) + ...
0.665 * v_sorted(25) +  0.243 * v_sorted(26) + ...
0.665 * v_sorted(27) +  0.243 * v_sorted(26) + ...
0.665 * v_sorted(27) +  0.243 * v_sorted(18));


V_avg_four_colom = 2 * (0.636 * v_sorted(21) +  0.243 * v_sorted(28) + ...
0.636 * v_sorted(23) +  0.243 * v_sorted(29) + ...
0.636 * v_sorted(25) +  0.243 * v_sorted(30));


V_avg_five_colom =2*( 0.665 * v_sorted(31) +  0.243 * v_sorted(28) + ...
    0.665 * v_sorted(31) +  0.243 * v_sorted(32)) + ...
    2*(0.665 * v_sorted(33) +  0.243 * v_sorted(32) + ...
    0.665 * v_sorted(33) +  0.243 * v_sorted(29) + ...
    0.665 * v_sorted(33) +  0.243 * v_sorted(34)) + ...
    2*(0.665 * v_sorted(35) +  0.243 * v_sorted(34) + ...
    0.665 * v_sorted(35) +  0.243 * v_sorted(30));


V_avg_six_colom = 2*(0.665 * v_sorted(33) +  0.243 * v_sorted(36)) ;


V_avg_seven_colom = 2*(0.665 * v_sorted(37) +  0.243 * v_sorted(36));


V_avg = V_corner * Asubc_36 / Across + V_edge_A  * Asubc_31_35 / Across + V_edge_B  * Asubc_26_30 / Across+ V_avg_three_colom * Asubc_1_25 / Across+ V_avg_four_colom * Asubc_1_25 / Across+ V_avg_five_colom * Asubc_1_25 / Across+ ...
    + V_avg_six_colom * Asubc_1_25 / Across + V_avg_seven_colom * Asubc_1_25 / Across;
Flow_rate = V_avg;
exp = load('flow_rate.mat');
Flow_rate_exp = exp.text_test(3,1);
fprintf('实验测定的流量: %.4f\n', Flow_rate_exp);
fprintf('计算得到的流量: %.4f\n', Flow_rate);



%Pranjape
V_corner_pran =  0.662 * v_sorted(1) +  0.243 * v_sorted(2) + ...
    0.662 * v_sorted(13) +  0.243 * v_sorted(12);
V_edge_A_pran = 0.662 * v_sorted(3) +  0.243 * v_sorted(2) + ...
    0.662 * v_sorted(3) +  0.243 * v_sorted(4) + ...
    0.662 * v_sorted(5) +  0.243 * v_sorted(4) + ...
    0.662 * v_sorted(5) +  0.243 * v_sorted(6) + ...
    0.662 * v_sorted(7) +  0.243 * v_sorted(6) + ...
    0.662 * v_sorted(7) +  0.243 * v_sorted(8) + ...
    0.662 * v_sorted(9) +  0.243 * v_sorted(8) + ...
    0.662 * v_sorted(9) +  0.243 * v_sorted(10) + ...
    0.662 * v_sorted(11) +  0.243 * v_sorted(10) + ...
    0.662 * v_sorted(11) +  0.243 * v_sorted(12);



V_edge_B_pran = 2 * (0.662 * v_sorted(3) +  0.243 * v_sorted(14) + ...
    0.662 * v_sorted(5) +  0.243 * v_sorted(15) + ...
    0.662 * v_sorted(7) +  0.243 * v_sorted(16) + ...
    0.662 * v_sorted(9) +  0.243 * v_sorted(17) + ...
    0.662 * v_sorted(11) +  0.243 * v_sorted(18));


V_avg_three_colom_pran = 2*( 0.662 * v_sorted(19) +  0.243 * v_sorted(14) + ...
0.662 * v_sorted(19) +  0.243 * v_sorted(20) + ...
0.662 * v_sorted(21) +  0.243 * v_sorted(20) + ...
0.662 * v_sorted(21) +  0.243 * v_sorted(15) + ...
0.662 * v_sorted(21) +  0.243 * v_sorted(22) + ...
0.662 * v_sorted(23) +  0.243 * v_sorted(22) + ...
0.662 * v_sorted(23) +  0.243 * v_sorted(16) + ...
0.662 * v_sorted(23) +  0.243 * v_sorted(24) + ...
0.662 * v_sorted(25) +  0.243 * v_sorted(24) + ...
0.662 * v_sorted(25) +  0.243 * v_sorted(17) + ...
0.662 * v_sorted(25) +  0.243 * v_sorted(26) + ...
0.662 * v_sorted(27) +  0.243 * v_sorted(26) + ...
0.662 * v_sorted(27) +  0.243 * v_sorted(18));


V_avg_four_colom_pran = 2 * (0.662 * v_sorted(21) +  0.243 * v_sorted(28) + ...
0.662 * v_sorted(23) +  0.243 * v_sorted(29) + ...
0.662 * v_sorted(25) +  0.243 * v_sorted(30));


V_avg_five_colom_pran =2*( 0.662 * v_sorted(31) +  0.243 * v_sorted(28) + ...
    0.662 * v_sorted(31) +  0.243 * v_sorted(32)) + ...
    2*(0.662 * v_sorted(33) +  0.243 * v_sorted(32) + ...
    0.662 * v_sorted(33) +  0.243 * v_sorted(29) + ...
    0.662 * v_sorted(33) +  0.243 * v_sorted(34)) + ...
    2*(0.662 * v_sorted(35) +  0.243 * v_sorted(34) + ...
    0.662 * v_sorted(35) +  0.243 * v_sorted(30));


V_avg_six_colom_pran = 2*(0.662 * v_sorted(33) +  0.243 * v_sorted(36)) ;


V_avg_seven_colom_pran = 2*(0.662 * v_sorted(37) +  0.243 * v_sorted(36));




V_avg_pran = V_corner_pran * Asubc_36 / Across + V_edge_A_pran  * Asubc_31_35 / Across + V_edge_B_pran  * Asubc_26_30 / Across+ V_avg_three_colom_pran * Asubc_1_25 / Across+ V_avg_four_colom_pran * Asubc_1_25 / Across+ V_avg_five_colom_pran * Asubc_1_25 / Across+ ...
    + V_avg_six_colom_pran * Asubc_1_25 / Across + V_avg_seven_colom_pran * Asubc_1_25 / Across;
Flow_rate_pran = V_avg_pran;
exp = load('flow_rate.mat');
Flow_rate_exp = exp.text_test(3,1);
fprintf('实验测定的流量: %.4f\n', Flow_rate_exp);
fprintf('计算得到的流量: %.4f\n', Flow_rate_pran);

% **图 1: 绘制散点图**
figure;
scatter(unique_points(:,1), unique_points(:,2), 50, v_avg, 'filled');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('完整的流场数据分布（去重后）');
grid on;
axis equal;

% =====================================
% ✅ 在图3点上补壁面点 + 遗传算法优化 m
% =====================================
X_exp = unique_points(:,1);
Y_exp = unique_points(:,2);
U_exp = v_avg;

% 通道尺寸
W_channel = 6 * P;
H_channel = 6 * P;

% 增加壁面点（速度为0）
n_wall = 20;
W_half = W_channel / 2;
H_half = H_channel / 2;

x_line = linspace(-W_half, W_half, n_wall)';
y_line = linspace(-H_half, H_half, n_wall)';

% 上下边界
Xw1 = [x_line; x_line];
Yw1 = [H_half*ones(n_wall,1); -H_half*ones(n_wall,1)];

% 左右边界
Xw2 = [W_half*ones(n_wall,1); -W_half*ones(n_wall,1)];
Yw2 = [y_line; y_line];

% 合并所有壁面点
X_wall = [Xw1; Xw2];
Y_wall = [Yw1; Yw2];
U_wall = zeros(size(X_wall));  % 速度设为0

% 拼接原始点 + 壁面点
X_aug = [X_exp; X_wall];
Y_aug = [Y_exp; Y_wall];
U_aug = [U_exp; U_wall];

% 归一化
r_x = abs(X_aug) / (W_channel/2);
r_y = abs(Y_aug) / (H_channel/2);
Umax = max(U_exp);

% 优化目标函数
error_fun = @(m) sum((Umax * (1 - r_x).^(1/m) .* (1 - r_y).^(1/m) - U_aug).^2);

% 遗传算法设置
opts = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 1000);
[m_opt, fval] = ga(error_fun, 1, [], [], [], [], 10, 40, [], opts);

% 输出结果
fprintf('\n🎯 遗传算法优化结果（含壁面点）：\n');
fprintf('最优幂律指数 m = %.6f\n', m_opt);
fprintf('拟合误差平方和 = %.6f\n', fval);

% 拟合速度
U_fit = Umax * (1 - r_x).^(1/m_opt) .* (1 - r_y).^(1/m_opt);

% 可视化拟合对比图
figure;
scatter3(X_aug, Y_aug, U_aug, 50, 'r', 'filled'); hold on;
scatter3(X_aug, Y_aug, U_fit, 50, 'b');
legend('实验 + 壁面点', '拟合速度');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('速度');
title(sprintf('含壁面点的幂律拟合结果：m = %.5f', m_opt));
grid on;
view(45, 30);


% ========================
% ✅ 计算相对误差指标
% ========================
relative_errors = abs((U_fit - U_aug) ./ U_aug);

% 避免除以0带来的Inf或NaN
relative_errors(U_aug == 0) = 0;

mean_rel_error = mean(relative_errors);
max_rel_error = max(relative_errors);

fprintf('\n 拟合误差评估：\n');
fprintf('平均相对误差（Mean Relative Error）= %.4f (%.2f%%)\n', mean_rel_error, mean_rel_error*100);
fprintf('最大相对误差（Max Relative Error） = %.4f (%.2f%%)\n', max_rel_error, max_rel_error*100);
% 计算相对误差
relative_errors = abs((U_fit - U_aug) ./ U_aug);
relative_errors(U_aug == 0) = 0;
mean_rel_error = mean(relative_errors);
max_rel_error = max(relative_errors);

% 显示图示（含误差信息）
title(sprintf('幂律拟合：m = %.5f | 平均相对误差 = %.2f%%', m_opt, mean_rel_error*100));


