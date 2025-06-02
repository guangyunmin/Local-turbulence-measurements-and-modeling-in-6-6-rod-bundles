
clc;
clear;
close all;

% =========================
% 🎯 Calculate Subchannel Area
% =========================
Umax = 1.23*1.087; % 你可以在这里输入任何你想要的 Umax 值或分布
D = 10;

Across = (100^2-(36*pi*(D/2)^2))/4;

% **读取并处理流场数据**
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

x = data(:,1);
y = data(:,2);
v = data(:,3);

% 生成对称点
x_full = [x; -x; x; -x; y; -y; y; -y];
y_full = [y; y; -y; -y; x; x; -x; -x];
v_full = [v; v; v; v; v; v; v; v];

% 去重
[unique_points, ~, idx] = unique([x_full, y_full], 'rows');
v_avg = accumarray(idx, v_full, [], @mean);

% 插值网格
[Xq, Yq] = meshgrid(linspace(-50, 50, 100), linspace(-50, 50, 100));
Vq = griddata(unique_points(:,1), unique_points(:,2), v_avg, Xq, Yq, 'cubic');
Z_min = min(Vq(:));


% 图1：原始散点
figure;
scatter(x, y, 50, v, 'filled');
colorbar;
xlabel('X (mm)'); ylabel('Y (mm)');
title('排序前的1/4流场');
for j = 1:length(x)
    text(x(j), y(j), num2str(v(j), '%.2f'), 'FontSize', 8, 'Color', 'k', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
grid on; view(-90, 90); colormap('jet');

% 图2：排序后
sortedData = sortrows([data(:,1), data(:,2), data(:,3)], [-2, -1]);
x_sorted = sortedData(:,1); y_sorted = sortedData(:,2); v_sorted = sortedData(:,3);
figure;
scatter(x_sorted, y_sorted, 50, v_sorted, 'filled');
colorbar;
xlabel('Y (mm)'); ylabel('X (mm)');
title('排序后的1/4流场');
grid on; view(-90, 90); colormap('jet');
for i = 1:length(x_sorted)
    text(x_sorted(i), y_sorted(i),sprintf('%.2f', v_sorted(i)), 'FontSize', 8, 'Color', 'k', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

% 图3：完整对称流场（实验点）
figure;
scatter(unique_points(:,1), unique_points(:,2), 50, v_avg, 'filled');
colorbar;
xlabel('X (mm)'); ylabel('Y (mm)');
title('完整的流场数据分布（去重后）');
grid on; axis equal; colormap('jet');

% ================== 计算混合长度 (lm) 的分布 ==================
C1 = 0.14;  
C2 = 0.08; 
C3 = 0.06;  
C4 = 26;
D = 0.01;  % Rod diameter (m)
P = 0.0167;% Rod bundle pitch (m)
R = D / 2; % Rod radius
m = 1/29.2;
Re_values= [7268.072,10884.969, 13055.358, 15961.487,20309.744];

Re_t_values = 1.35 .* 0.4 .* Re_values./log(Re_values)/2;
decay_factor = Re_t_values(5) / C4;
channel_width = 0.1;
channel_height = 0.1;

% Set up computational grid (-3P, 3P) to ensure all rods are included
Nx = 300; 
Ny = 300;
x = linspace(-3*P, 3*P, Nx);
y = linspace(-3*P, 3*P, Ny);
[X, Y] = meshgrid(x, y);

% Calculate rod center positions, aligning the center to (0,0)
num_rows = 6; % 6x6 rod bundle
num_cols = 6;
dx = P; 
dy = P;

x_centers = ((1:num_cols) - (num_cols + 1) / 2) * dx;
y_centers = ((1:num_rows) - (num_rows + 1) / 2) * dy;
[Xc, Yc] = meshgrid(x_centers, y_centers);
rod_centers = [Xc(:), Yc(:)];

% Compute the distance to the nearest rod center
min_distances = inf(size(X));
for k = 1:size(rod_centers, 1)
    distances = sqrt((X - rod_centers(k,1)).^2 + (Y - rod_centers(k,2)).^2);
    min_distances = min(min_distances, distances);
end

% Compute the mixing length (lm), setting regions with lm < R to NaN
lm = zeros(size(X));
idx_core = (min_distances >= R);
lm(idx_core) = (C1 - C2 * (1 - (-R + min_distances(idx_core)) / (max(min_distances(:)) - R)).^2 ...
                - C3 * (1 - (-R + min_distances(idx_core)) / (max(min_distances(:)) - R)).^4) ...
                .* (1 - exp(decay_factor * ((R - min_distances(idx_core)) / (max(min_distances(:)) - R))));

            
                     
% 可视化混合长度 lm 分布
figure;
surf(X, Y, lm,  'EdgeColor', 'none'); 
colorbar;
colormap(jet);
xlabel('\it{X}, -', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('\it{Y}, -', 'FontName', 'Times New Roman', 'FontSize', 12);
zlabel('\lambda_{rod}, -', 'Interpreter', 'tex', 'FontName', 'Times New Roman', 'FontSize', 12);
% 设置坐标轴属性
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 3, ...               % 坐标轴刻度线加粗
    'GridLineStyle', '-', ...
    'GridAlpha', 0.4, ...
    'GridColor', [0.3 0.3 0.3], ...
    'LineWidth', 1.5);                  % 网格线加粗
grid on;
xticks(-1:0.4:1);
yticks(-1:0.4:1);



% ================== 训练 RBF 神经网络 ==================
% 训练 RBF 神经网络，允许 Umax 作为输入变量
if ~exist('net', 'var')
    y_values = linspace(0, 1, 1000)'; 
    
    % 动态变化的 Umax (在训练阶段增加 Umax 变化的多样性)
    Umax_values = 0.3 + 2.5 * rand(length(y_values), 1); 
    
    % 计算混合长度 (lm) 与 u 的关系
    lm_function = @(y) 0.14 - 0.08 * (1 - y).^2 - 0.06 * (1 - y).^4;
    lm_values = lm_function(y_values);
    
    % 目标函数 u = Umax * y^(1/7)
    u_function = @(y, Umax) Umax .* y.^(1/(6.32));
    u_values = u_function(y_values, Umax_values);

    % 训练 RBF 神经网络
    input_data = [lm_values, Umax_values]';
    u_train = u_values';

    goal = 1e-6;  
    spread = 0.1; 
    net = newrb(input_data, u_train, goal, spread, 100);
end

% ================== 预测 u 分布 (使用指定的 Umax) ==================
% 这里直接指定 Umax，例如设为 1.2（你可以自行修改）


% 构建预测输入
lm_values = lm(:);
Umax_values = Umax * ones(size(lm_values)); % 使用指定的 Umax

% 预测输入与预测
prediction_input = [lm_values, Umax_values]';
u_pred = net(prediction_input);
u_pred = reshape(u_pred, size(lm));
% ================== 应用整个通道的幂律速度分布调制 ==================

% 通道尺寸和幂律参数
W = 6 * P;  % 这里选择流道宽度为 6P，与混合长度网格对齐
H = 6 * P;  % 同理

% 为了配合混合长度计算的 X, Y 网格，重新归一化坐标
r_x = abs(X) / (W/2);
r_y = abs(Y) / (H/2);

% 通道幂律分布（使用均匀最大值归一化）
U_channel_x = (1 - r_x).^m;
U_channel_y = (1 - r_y).^m;
U_channel = Umax.*U_channel_x .* U_channel_y;
U_channel(U_channel < 0) = 0;
U_original = U_channel;
% ========= 面积平均速度（缩放前） =========
dx = W / (Nx - 1);
dy = H / (Ny - 1);
dA = dx * dy;
total_area = W * H;
U_avg_original = sum(U_original(:) * dA) / total_area;
% ========= 缩放以归一化平均速度为 1 =========
scale_factor = 1 / U_avg_original;
U_scaled = U_original * scale_factor;
% ========= 缩放后再计算一次平均值以验证 =========
U_avg_scaled = sum(U_scaled(:) * dA) / total_area;
% ---------- 可视化原始分布 ----------
figure;
surf(X, Y, U_scaled, 'EdgeColor', 'none');
idx=find(Y>=0);
XXX=X(idx)*1000;
YYY=Y(idx)*1000;
ZZZ=U_original(idx);
XXX(:);
YYY(:);
ZZZ(:);







colormap(jet); colorbar;
title('原始幂律速度分布 (n = 1)');
xlabel('X'); ylabel('Y'); zlabel('Velocity');
shading interp;
% 将原始预测速度乘以幂律通道修正项
u_pred = u_pred .* U_scaled;

% ================== 可视化预测结果 ==================

% 2. 预测 u 的 3D 表面图
figure;
surf(X, Y, u_pred, 'EdgeColor', 'none');
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Predicted u');
title(['3D Predicted u Distribution with Specified Umax = ', num2str(Umax)]);
view(3);
grid on;


% 计算面积加权平均流速
x_centers = (x(1:end-1) + x(2:end)) / 2;
y_centers = (y(1:end-1) + y(2:end)) / 2;
[X_centers, Y_centers] = meshgrid(x_centers, y_centers);

% 插值计算网格中心的流速 (使用 u_pred 作为速度场)
v_centers = interp2(X, Y, u_pred, X_centers, Y_centers, 'linear');
v_centers(isnan(v_centers)) = 0;

% 计算网格中心的面积
dx = x_centers(2) - x_centers(1);
dy = y_centers(2) - y_centers(1);
cell_area = dx * dy;

% 计算网格中心到最近杆中心的距离 (用于判断流体区域)
min_distances_centers = inf(size(X_centers));
for k = 1:size(rod_centers, 1)
    distances_centers = sqrt((X_centers - rod_centers(k,1)).^2 + (Y_centers - rod_centers(k,2)).^2);
    min_distances_centers = min(min_distances_centers, distances_centers);
end

% 生成流体区域掩膜 (距离大于杆半径 R 的区域为流体区域)
fluid_mask = min_distances_centers >= R;

% 计算面积加权平均流速
area_weighted_sum = sum(v_centers(fluid_mask) * cell_area);
total_area = sum(fluid_mask(:)) * cell_area;
area_average = area_weighted_sum / total_area;

% 计算最大流速 (仅在流体区域内)
max_velocity = max(v_centers(fluid_mask));

% 输出计算结果
disp(['面积加权平均流速: ', num2str(area_average)]);
disp(['最大流速: ', num2str(max_velocity)]);


% =================== 添加平滑过渡处理的代码 ===================
% 设置固体区域流场为0，并且保证平滑过渡
transition_width = 1e-3; % 过渡区域的宽度 (可根据需要调整)

% 计算到最近固体区域的距离
distance_to_solid = min_distances_centers - R;

% 生成平滑过渡系数 (0 表示在固体区域，1 表示完全在流体区域)
smooth_transition = ones(size(distance_to_solid));
smooth_transition(distance_to_solid < 0) = 0; % 完全在固体区域为0

% 在过渡区域内进行平滑过渡 (采用线性插值)
transition_idx = (distance_to_solid >= 0) & (distance_to_solid <= transition_width);
smooth_transition(transition_idx) = ...
    distance_to_solid(transition_idx) / transition_width;

% 插值 smooth_transition 到与 u_pred 相同的网格大小 (保持所有区域有值)
smooth_transition_resized = interp2(X_centers, Y_centers, smooth_transition, X, Y, 'linear', 0);

% 应用平滑过渡到预测的速度场 u_pred
u_pred = u_pred .* smooth_transition_resized;

% 对速度场 (U, V) 也应用平滑过渡
[U, V] = gradient(u_pred, x, y);
U = U .* smooth_transition_resized;
V = V .* smooth_transition_resized;


% 确保 u_pred 中固体区域为 0 而不是 NaN
u_pred(isnan(u_pred)) = 0;

% 对速度场 (U, V) 也应用相同的处理，避免在 quiver 中出现空白
U(isnan(U)) = 0;
V(isnan(V)) = 0;
% ================== 可视化更新后的流场 ==================
figure;

hold on;
surf(X, Y, u_pred,  'EdgeColor', 'k'); 
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');
title(['Flow Field with Smooth Transition (Umax = ', num2str(Umax), ')']);
view(3)

% 找到 u_pred 中的最大值及其坐标
[max_value, max_index] = max(u_pred(:));
[max_i, max_j] = ind2sub(size(u_pred), max_index);
max_x = X(max_i, max_j);
max_y = Y(max_i, max_j);

% 在图中标出最大值的位置，并显示数值
plot3(max_x, max_y,max_value, 'k*', 'MarkerSize', 10, 'LineWidth', 2); % 红色星形标记最大值
text(max_x, max_y,max_value+0.1, sprintf('Max: %.2f', max_value), 'Color', 'w', ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
hold off;

XX=X(:);
YY=Y(:);
ZZ=u_pred(:);

% 找出 y = 0 上的预测流速
target_y = 0 / 1000;           % 目标 Y 值
tolerance = 0.0005;                 % 容差
y_idx = find(abs(Y(:,1) - target_y) < tolerance);  % 在第一列中查找
u_pred_y0 = u_pred(y_idx, :);         % 提取 y = 0 对应行的预测速度
x_line = x;                                % x 方向坐标

% 可视化 y = 0 上的流速分布
figure;
plot(x_line, u_pred_y0, 'b-', 'LineWidth', 2);
xlabel('X (m)');
ylabel('u at Y = 0');
title('Predicted u Distribution along Y = 0');
grid on;

% ================= 误差分析（基于图3实验点） =================
x_exp = unique_points(:,1) / 1000;  % mm -> m
y_exp = unique_points(:,2) / 1000;
v_exp = v_avg;

u_interp = interp2(X, Y, u_pred, x_exp, y_exp, 'linear');
valid = ~isnan(u_interp);
v_exp = v_exp(valid);
u_interp = u_interp(valid);

relative_error = abs(u_interp - v_exp) ./ abs(v_exp);
mean_error = mean(relative_error);

fprintf('图3实验点的平均相对误差为：%.4f (%.2f%%)\n', mean_error, mean_error * 100);

% ================== 可视化实验值与预测值的 3D 对比图 ==================
figure;
scatter3(x_exp * 1000, y_exp * 1000, v_exp, 50, 'r', 'filled'); hold on;
scatter3(x_exp * 1000, y_exp * 1000, u_interp, 50, 'b');  % 坐标单位转为 mm
legend('实验速度', '预测速度');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('速度');
title('实验点与预测速度的3D空间分布对比');
grid on;
view(45, 30);

% ================== 可视化更新后的流场 ==================
figure;
surf(X, Y, u_pred); 
ZZ=(u_pred(:));
hold on;
scatter3(x_exp, y_exp, v_exp, 60, 'g', 'filled');
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');



