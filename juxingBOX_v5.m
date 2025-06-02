clc; clear; close all;

% **读取并处理流场数据**
% 你的数据点
data = [
0	0	1.3157788744721
0	8.3333	1.28150498910915
0	25	1.21517117223561
0	16.6666666666667	1.27225753857891
-8.33333333333333	16.6666666666667	1.20737254930127
-16.6666666666667	16.6666666666667	1.21741400728532
-16.6667	25	1.19363302832612
8.33333333333333	16.6666666666667	1.18724833863244
16.6666666666667	16.6666666666667	1.19981239385144
16.6667	25	1.19305975903111
16.6667	41.6667	1.14269563124375
0	33.3333333333333	1.25527811331408
-8.33333333333333	33.3333333333333	1.19706370329333
-16.6666666666667	33.3333333333333	1.22026377974166
-25	33.3333333333333	1.18409995709111
-33.3333333333333	33.3333333333333	1.19462538680417
8.33333333333333	33.3333333333333	1.17744340279046
16.6666666666667	33.3333333333333	1.19660715160325
25	33.3333333333333	1.19082893136763
33.3333333333333	33.3333333333333	1.17398172253408
33.3333	41.6667	1.14496847346418
0	41.6667	1.17399165351035
-16.6667	41.6667	1.1407653732544
-33.3333	41.6667	1.1317058433363
0	48.335	1.20212446073426
-8.33333333333333	48.335	1.11777303308902
-16.6666666666667	48.335	1.17990097101792
-25	48.335	1.13861006160721
-33.3333333333333	48.335	1.17825137965713
-41.6666666666667	48.335	1.08384054703005
8.33333333333333	48.335	1.12444891107579
16.6666666666667	48.335	1.16417402213646
25	48.335	1.14538593865467
33.3333333333333	48.335	1.1599233083827
41.6666666666667	48.335	1.07858535719249
-48.335	48.335	1.10378878645276
48.335	48.335	1.11895367911525];
save('data5.mat','data')
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
fprintf('Cc_1_25: %.4f\n', Cc_1_25);
fprintf('Ce_1_25: %.4f\n', Ce_1_25);

fprintf('Cc_26_30: %.4f\n', Cc_26_30);
fprintf('Ce_26_30: %.4f\n', Ce_26_30);

fprintf('Cc_31_35: %.4f\n', Cc_31_35);
fprintf('Ce_31_35: %.4f\n',Ce_31_35);

fprintf('Cc_36: %.4f\n', Cc_36);
fprintf('Ce_36: %.4f\n', Ce_36);


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
Flow_rate_exp = exp.text_test(5,1);
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
Flow_rate_exp = exp.text_test(5,1);
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

% **图 2: 流场云图 + 透明实体圆柱**
figure;
hold on;
surf(Xq, Yq, Vq, 'EdgeColor', 'none'); % 画出流场 (不透明)
shading interp;
colorbar;
Xq = Xq(:);
Yq = Yq(:);
Vq = Vq(:);

valid_idx = ~isnan(Vq);  
Xq = Xq(valid_idx);
Yq = Yq(valid_idx);
Vq = Vq(valid_idx);


xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Velocity');
title('流场云图 + 透明圆柱阵列');
grid on;
view(3);

% **在 `surf` 图上叠加透明实体圆柱**
% **设置圆柱参数**
D = 10;   % 圆柱直径 (mm)
P = 16.7; % 圆柱阵列间距 (mm)
R = D / 2; % 圆柱半径

% **计算圆柱中心位置**
num_rows = 6;
num_cols = 6;
dx = P; 
dy = P;

x_centers = ((1:num_cols) - (num_cols + 1) / 2) * dx;
y_centers = ((1:num_rows) - (num_rows + 1) / 2) * dy;
[Xc, Yc] = meshgrid(x_centers, y_centers);
rod_centers = [Xc(:), Yc(:)];

AA=Xc(:);
BB=Yc(:);
% **绘制透明实体圆柱**
theta = linspace(0, 2*pi, 50); % 圆的角度
Z_bottom = Z_min; % 圆柱从流场 Z 的最小值开始
Z_top = Z_min + 0.15; % 圆柱高度 0.3
for i = 1:size(rod_centers,1)
    % 生成圆柱的坐标
    [Xcyl, Ycyl, Zcyl] = cylinder(R, 50);
    Xcyl = Xcyl + rod_centers(i,1); % 移动到正确的位置
    Ycyl = Ycyl + rod_centers(i,2);
    Zcyl = Zcyl * (Z_top - Z_bottom) + Z_bottom; % 设定圆柱高度

    % 画出透明实体圆柱
    surf(Xcyl, Ycyl, Zcyl, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'k'); 
end
colormap('jet')
hold off;
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

