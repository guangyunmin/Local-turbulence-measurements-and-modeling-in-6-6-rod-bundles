clc; clear; close all;

% ================== 计算混合长度 (lm) 的分布 ==================
C5 =3.81;
% Umax = 1.154*0.389; % 你可以在这里输入任何你想要的 Umax 值或分布
Umax = 1.2412*0.389; % 你可以在这里输入任何你想要的 Umax 值或分布
% Set parameters
C1 = 0.14;  
C2 = 0.08; 
C3 = 0.06;
C4 = 26;
D = 0.01;  % Rod diameter (m)
P = 0.0167;% Rod bundle pitch (m)
R = D / 2; % Rod radius

Re_values= [7268.072,10884.969, 13055.358, 15961.487,20309.744];

Re_t_values = 1.35 .* 0.4 .* Re_values./log(Re_values);
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
xlabel('X (m)');
ylabel('Y (m)');
title('Mixing Length (lm) Distribution (lm < R regions as blank)');

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
    u_function = @(y, Umax) Umax .* y.^(1/(6.535));
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



% Define known y positions (in meters) and corresponding coefficients
y_known = [0.0, 8.3, 16.7, 25.0, 33.3, 41.7, 48.3] * 1e-3;

% Coefficients for v(x) = a0 + a1 * cos(f1 * x) + ... 
coefficients = [
    0.5184, 0.0174, 0.0071, 0.0004, -0.0082, -0.0008, 0.0025;
    0.4998, 0.0243, 0.0044, 0.0107,  0.0002,  0.0022, -0.0010;
    0.4966, 0.0120, -0.0015, 0.0108, 0.0012, -0.0009,  0.0090;
    0.4965, -0.0001, -0.0019, 0.0050, -0.0016, 0.0010, -0.0007;
    0.4998, 0.0083, -0.0064, 0.0025, 0.0018,  0.0024,  0.0066;
    0.4766, 0.0123,  0.0027, 0.0137, -0.0011, 0.0025, -0.0014;
    0.4769, 0.0121, -0.0144, 0.0077,  0.0020,  0.0042,  0.0121;
];




% Cosine frequencies
cos_frequencies = [0.0650, 0.1300, 0.1950, 0.2600, 0.3250, 0.3900];

% Initialize v_decay_matrix
v_decay_matrix = zeros(size(X));

% 双向插值，保证 X 和 Y 方向都平滑衰减
for j = 1:Nx
    for i = 1:Ny
        y_pos = abs(Y(i, j));
        x_pos = abs(X(i, j));
        
        % 分别对 X 和 Y 方向进行插值，保证双向对称性
        interpolated_coeffs_x = arrayfun(@(k) interp1(y_known, coefficients(:, k), x_pos, 'spline'), 1:size(coefficients, 2));
        interpolated_coeffs_y = arrayfun(@(k) interp1(y_known, coefficients(:, k), y_pos, 'spline'), 1:size(coefficients, 2));
        
        v_decay_x = interpolated_coeffs_x(1);
        v_decay_y = interpolated_coeffs_y(1);
        
        for k = 1:length(cos_frequencies)
            v_decay_x = v_decay_x + interpolated_coeffs_x(k+1) * cos(cos_frequencies(k) * x_pos);
            v_decay_y = v_decay_y + interpolated_coeffs_y(k+1) * cos(cos_frequencies(k) * y_pos);
        end
        
        % 乘积模式保证双向相同的衰减规律
        v_decay_matrix(i, j) = v_decay_x * v_decay_y;
    end
end

% Normalize v_decay and apply it to lm
% v_decay = (v_decay - min(v_decay(:))) / (max(v_decay(:)) - min(v_decay(:)));
% scale_min = 0.97;
% scale_max = 1.0;
% v_decay = scale_min + (scale_max - scale_min) * ...
%     ( v_decay_matrix - min( v_decay_matrix(:))) / (max( v_decay_matrix(:)) - min( v_decay_matrix(:)));
% 
u_pred = C5*(u_pred .* v_decay_matrix);
% u_pred = u_pred .* v_decay;

% 在 X 和 Y 方向的 -50 到 -48.335 以及 48.335 到 50 的区域中平滑过渡到0
transition_start = 48.335e-3;
transition_end = max(X(:));

% 对 X 方向的平滑过渡处理
for j = 1:Nx
    x_pos = abs(x(j));
    if x_pos >= transition_start && x_pos <= transition_end
        % 使用线性衰减，保证50处为0，48.335处平滑
        transition_factor = (transition_end - x_pos) / (transition_end - transition_start);
        u_pred(:, j) = u_pred(:, j) * transition_factor;
    end
end

% 对 Y 方向的平滑过渡处理
for i = 1:Ny
    y_pos = abs(y(i));
    if y_pos >= transition_start && y_pos <= transition_end
        % 使用线性衰减，保证50处为0，48.335处平滑
        transition_factor = (transition_end - y_pos) / (transition_end - transition_start);
        u_pred(i, :) = u_pred(i, :) * transition_factor;
    end
end 




% ================== 可视化预测结果 ==================
% 1. 预测 u 分布的等高线图
figure;
contourf(X, Y, u_pred, 50, 'LineColor', 'none'); 
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');
title(['Predicted u Distribution with Specified Umax = ', num2str(Umax)]);
axis equal;

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

% 3. 预测流场图 (基于 u 分布的速度场)
[U, V] = gradient(u_pred, x, y); 
figure;
quiver(X, Y, U, V, 'k'); 
hold on;
contourf(X, Y, u_pred, 50, 'LineColor', 'none'); 
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');
title(['Flow Field with Specified Umax = ', num2str(Umax)]);
axis equal;

 



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




figure;
contourf(X, Y, u_pred, 50, 'LineColor', 'none'); 
ZZ=u_pred(:);
colorbar;
colormap(jet);
xlabel('X (m)');
ylabel('Y (m)');
title(['Predicted u Distribution with Specified Umax = ', num2str(Umax)]);
axis equal;
