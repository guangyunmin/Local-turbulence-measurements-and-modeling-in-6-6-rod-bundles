clc; clear; close all; 

% **Parameter Definition**
P = 16.7;     % Structure size parameter
D = 10;       % Diameter of the circular region
f_g = 0.98;   % Base flow coefficient
f_0 = 1;   % Maximum flow coefficient
m = 7;        % Power exponent
W = 1.665;  
R0 = D /2;
Rnew = R0;   % Initial radius


% **Create polar coordinate grid**
theta = linspace(0, pi / 4, 50);  
r = linspace(0, ((P/2-W) * (sqrt(2))), 50); 
[THETA, R] = meshgrid(theta, r);

% **Compute flow field**
sec_theta = 1 ./ cos(THETA);
f_c = f_g + (f_0 - f_g) .* (1 - ((sqrt(2) - sec_theta) / (sqrt(2) - 1)).^m);
denominator = ((P / 2 - W) .* sec_theta - Rnew);
F = f_c .* (1 - (1 - (R - R0) ./ denominator).^m);


F(R <= Rnew) = 0; 

% **Convert to Cartesian coordinates**
[X, Y] = pol2cart(THETA, R);
X_data = [X(:), Y(:)];
F_data = F(:);



% **Compute P2 structured grid**
theta_hou = atan((P/2 - W) / (P/2 - W));
arc_theta_hou = linspace(theta_hou, 0, 30);
arc_theta_qian = linspace(pi/4, 0, 30); 
R_qian = Rnew; 
R_hou = Rnew; 
arc_x_qian = R_qian * cos(arc_theta_qian);  
arc_y_qian = R_qian * sin(arc_theta_qian);
arc_x_hou = R_hou * cos(arc_theta_hou);
arc_y_hou = R_hou * sin(arc_theta_hou);
arc_points_hou = [arc_x_hou', arc_y_hou'];
arc_points_qian = [arc_x_qian', arc_y_qian'];
B = [P/2 - W, 0];                  
C = [P/2 - W, P/2 - W];                
D_hou = [R_hou * cos(theta_hou) , R_hou * sin(theta_hou)];
D_qian = [R_qian * cos(pi/4) , R_qian * sin(pi/4)];
P1_points = [arc_points_qian; B; C; D_qian];
P2_points = [arc_points_hou; B; C; D_hou];

% **Plot P1 and P2**
figure(1); 
hold on; 
axis equal; 
grid on;
% **Plot P1 (Green)**
plot(P1_points(:,1), P1_points(:,2), 'g-', 'LineWidth', 1.5);
scatter(P1_points(:,1), P1_points(:,2), 10, 'g', 'filled');

% **Plot P2 (Red)**
plot(P2_points(:,1), P2_points(:,2), 'r-', 'LineWidth', 1.5);
scatter(P2_points(:,1), P2_points(:,2), 10, 'r', 'filled');

% **Label key points**
text(B(1), B(2), ' B', 'FontSize', 12, 'Color', 'black');
text(C(1), C(2), ' C', 'FontSize', 12, 'Color', 'black');
text(D_qian(1), D_qian(2), ' D_qian', 'FontSize', 12, 'Color', 'black');
text(D_hou(1), D_hou(2), ' D_hou', 'FontSize', 12, 'Color', 'black');

xlabel('X Axis');
ylabel('Y Axis');
title('P1 (Green) and P2 (Red) Regions');
legend('P1 Boundary', 'P1 Points', 'P2 Boundary', 'P2 Points');
hold off;

P1_area = polyarea(P1_points(:,1), P1_points(:,2));
P2_area = polyarea(P2_points(:,1), P2_points(:,2));
juxing_area = W * (P / 2 - W);

% **Output the areas of P1 and P2**
fprintf('P1 structured grid area: %.4f\n', P1_area);
fprintf('P2 structured grid area: %.4f\n', P2_area);
fprintf('Rectangle area: %.4f\n', juxing_area);
fprintf('juxing_area + P2_area: %.4f\n', juxing_area + P2_area);


% **Create P2 structured grid**
P2_x = linspace(min(P2_points(:,1)), max(P2_points(:,1)), 300);
P2_y = linspace(min(P2_points(:,2)), max(P2_points(:,2)), 300);
[P2_X, P2_Y] = meshgrid(P2_x, P2_y);

[in, ~] = inpolygon(P2_X, P2_Y, P2_points(:,1), P2_points(:,2));
P2_X(~in) = NaN;
P2_Y(~in) = NaN;

% % **Mask out values where X > (P/2 - W)**
% mask = (X <= (P / 2 - W));
% X(~mask) = NaN;
% Y(~mask) = NaN;
% F(~mask) = NaN;


% **RBF interpolation**
X_train = [X(:), Y(:)];
F_train =  F_data;
valid_idx = ~isnan(X_train(:,1)) & ~isnan(X_train(:,2));
X_train = X_train(valid_idx, :);
F_train = F_train(valid_idx);

X_target = [P2_X(:), P2_Y(:)];
valid_target_idx = ~isnan(X_target(:,1)) & ~isnan(X_target(:,2));
X_target = X_target(valid_target_idx, :);

D_train = pdist2(X_train, X_train);
sigma = 0.16*median(D_train(:));  
rbf_func = @(r) exp(-(r / sigma).^2);
Phi_train = rbf_func(D_train);
lambda = 1e-6;  
Phi_train = Phi_train + lambda * eye(size(Phi_train));
weights = Phi_train \ F_train;

D_target = pdist2(X_target, X_train);
Phi_target = rbf_func(D_target);

F_P2_grid = NaN(size(P2_X));
F_P2_grid(valid_target_idx) = Phi_target * weights;

% **Right-side rectangular flow field `F_R`**
Nx_R = 300; Ny_R = 300;
x_R = linspace(P/2 - W, P/2, Nx_R);  
y_R = linspace(0, P/2 - W, Ny_R);        
[X_R, Y_R] = meshgrid(x_R, y_R);

% **Ensure valid values at P/2-W in `F_P2_grid`**
[~, idx] = min(abs(P2_X(1,:) - P/2 -W)); 
F_P2_P_boundary = F_P2_grid(:, idx);  % Directly extract values at X=P/2-W
F_P2_P_boundary = fillmissing(F_P2_P_boundary, 'linear'); % Linearly fill NaN values

% **Compute `F_R`, ensuring continuity at X = P/2**
F_R = NaN(size(X_R));

for i = 1:length(x_R)
    for j = 1:length(y_R)
        if X_R(j, i) == P/2-W
            F_R(j, i) = F_P2_P_boundary(j); % Ensure continuity at P/2
        else
            % **Extract boundary values from left P2 region**
            Fy = F_P2_P_boundary(j);
            
            % **Apply power-law decay in the X direction**
            Fx = (1 - ((X_R(j, i) - (P/2-W)) / W).^7);
            
            % **Final flow field computation**
            F_R(j, i) = Fy .* Fx;
        end
    end
end

% **Merge P2 grid and right rectangular grid**
X_combined = [P2_X, X_R]; 
Y_combined = [P2_Y, Y_R]; 
F_combined = [F_P2_grid, F_R];

F_combined(isnan(F_combined)) = 0; % Avoid NaN affecting visualization

% **Figure 2: Flow field before interpolation (Cartesian coordinates)**
figure(2);
surf(X, Y, F, 'EdgeColor', 'none');
colormap('jet'); 
colorbar;
xlabel('X'); ylabel('Y'); zlabel('F');
title('Flow field before interpolation (Cartesian coordinates)');
xlim([0,P/2-W])
zlim([0,1])
caxis([0 1])
% **Figure 3: Interpolated complete flow field**
figure(3);
surf(P2_X, P2_Y, F_P2_grid, 'EdgeColor', 'none');
colormap('jet'); 
colorbar;
xlabel('X'); ylabel('Y'); zlabel('F');
title('Interpolated flow field');

figure(4);
hold on;

% 绘制散点图
consurf(X_data(:,1), X_data(:,2), 10, 'b', 'filled');  
scatter(P2_X(:), P2_Y(:), 10, 'r', 'filled');  
PXX=P2_X(:);
PYY=P2_Y(:);

% 将 NaN 替换为 0
PXX(isnan(PXX)) = 0;
PYY(isnan(PYY)) = 0;



% 设置轴标签
xlabel('\it X\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('\it Y\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');

% 设置 legend 并调整符号大小
h = legend('Background Grid', 'Foreground Grid', 'Location', 'northwest');
set(h, 'Box', 'off');                
set(h, 'FontName', 'Times New Roman'); 
set(h, 'FontSize', 14); 
set(gcf, 'Renderer', 'painters'); % 使用矢量渲染
h.ItemTokenSize = [20, 10]; % **放大 legend 符号**

daspect([1 1 0.5]); % X and Y axes equal, Z scaled dynamically

% **设置坐标轴格式，确保所有刻度可见**
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
box on; % **确保所有框线显示**

% **确保所有刻度都显示**
xticks(min(X_data(:,1)):1:max(X_data(:,1))); % X 轴刻度
yticks(min(X_data(:,2)):1:max(X_data(:,2))); % Y 轴刻度

% **设置轴颜色确保所有刻度都可见**
set(gca, 'XColor', 'k', 'YColor', 'k'); % 确保 X、Y 轴的刻度线显示
set(gca, 'TickDir', 'in');  % 让刻度向外

% 设置 x 和 y 轴范围
xlim([0 max(X_data(:,1))]);
ylim([0 max(X_data(:,2))]);

hold off;







% **Figure 5: Right-side flow field**
figure(5);
surf(X_R, Y_R, F_R, 'EdgeColor', 'none'); 
colormap(jet);
colorbar;
xlabel('X'); ylabel('Y'); zlabel('F_R');
title('Right-side flow field (Rectangular grid)');
grid on;
view(45, 30);

% **Figure 6: Merged flow field**(juxing_area + P2_area)
figure(6);
surf(X_combined, Y_combined, F_combined, 'EdgeColor', 'none'); 
colormap(jet);
colorbar;
xlabel('X'); ylabel('Y'); zlabel('Flow Field');
title('Merged flow field');
grid on;
view(30, 40);


% **Define triangular grid region (F_U)**
Nx_U = 300; 
Ny_U = 300;
y_U = linspace(P/2 - W, P/2, Ny_U);  % Upper region (P/2-W < Y <= P/2)
x_U = linspace(P/2 - W, P/2, Nx_U);  % Same X range as F_R

[X_U, Y_U] = meshgrid(x_U, y_U);

% **Mask out points where Y > X to keep only the triangular region**
mask_triangular = (Y_U <= X_U); 

% **Ensure valid values at Y = P/2-W from `F_R`**
[~, idx_U] = min(abs(Y_R(:,1) - (P/2 - W))); 
F_R_Y_boundary = F_R(idx_U, :);  % Extract values at Y = P/2-W
F_R_Y_boundary = fillmissing(F_R_Y_boundary, 'linear'); % Fill NaNs

% **Compute F_U using power-law decay in Y direction**
F_U = NaN(size(X_U));

for i = 1:length(x_U)
    for j = 1:length(y_U)
        if mask_triangular(j, i) % Only compute where Y <= X
            if Y_U(j, i) == P/2 - W
                F_U(j, i) = F_R_Y_boundary(i); % Maintain continuity
            else
                Fy = (1 - ((Y_U(j, i) - (P/2 - W)) / W).^m);
                F_U(j, i) = F_R_Y_boundary(i) .* Fy;
            end
        end
    end
end

% **Extract only the valid triangular region data**
X_tri = X_U(mask_triangular);
Y_tri = Y_U(mask_triangular);
F_tri = F_U(mask_triangular);

% **Merge the structured flow fields**juxing_area + P2_area + triangular_area
X_final = [X_combined(:); X_tri(:)]; 
Y_final = [Y_combined(:); Y_tri(:)]; 
F_final = [F_combined(:); F_tri(:)];

F_final(isnan(F_final)) = 0; % Avoid NaN affecting visualization
X_final(isnan(X_final)) = 0; % Avoid NaN affecting visualization
Y_final(isnan(Y_final)) = 0; % Avoid NaN affecting visualization


% **Figure 7: Final Flow Field
figure(7);
scatter3(X_final, Y_final, F_final, 10, F_final, 'filled'); 
colormap(jet);
colorbar;
% **Set labels with Times New Roman Italic**
xlabel('\it X\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('\it Y\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
zlabel('\it Z\rm, -', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
title('1 / 2 corner subchannel', 'FontName', 'Times New Roman', 'FontSize', 14);

% **Ensure X and Y axes have equal scaling but leave Z free**
daspect([1 1 0.5]); % X and Y axes equal, Z scaled dynamically

% **Set tick intervals to 1**
xticks(min(X_final):1:max(X_final));
yticks(min(Y_final):1:max(Y_final));
% % zticks(min( F_final):1:max( F_final));
xlim([0 max(X_final)]);
ylim([0 max(Y_final)]);
% zlim([0 max(F_final)]);
% **Set grid and view**
grid on;
view(30, 40);
% **设置 X, Y, Z 轴的线条更粗**
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
set(gcf, 'Renderer', 'painters'); % 使用矢量渲染




% **Find the maximum value of F_final and its index**
[max_F, max_idx] = max(F_final(:)); % Find the maximum value
max_X = X_final(max_idx); % Get the corresponding X coordinate
max_Y = Y_final(max_idx); % Get the corresponding Y coordinate
hold on;
plot3(max_X, max_Y, max_F, 'bo', 'MarkerSize', 14, 'MarkerFaceColor', 'b'); % Mark maximum value
text(max_X-0.4, max_Y+0.2, max_F+0.4,  sprintf('\\it f\\rm_{c}: %.2f', max_F), ...
    'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'r', 'Interpreter', 'tex');
scatter3(D/2+W, 0, 0.98, 200, [0.5, 1, 1], 'filled'); 
text(D/2+W, -2.4, 0.98+0.4,sprintf('\\it f\\rm_{g}: %.2f', 0.98), ...
    'FontSize', 14, 'FontName', 'Times New Roman', ...
    'Interpreter', 'tex', 'Color', 'r');

% **Add a cylinder (fuel rod) centered at (0,0)**


% **定义圆柱参数**
x_center = 0;  % 圆柱中心坐标 (0,0)
y_center = 0;
z_min = min(F_final(:)); % 圆柱底部
z_max = max(F_final(:)); % 圆柱顶部
radius = 5;  % 圆柱半径
num_points = 50;  % 控制圆柱平滑度
cyl_color = [0.8, 0.7, 1]; % **浅紫色**
alpha_value = 0.5; % **透明度 (0=完全透明, 1=完全不透明)**

% **创建圆柱的外表面**
theta = linspace(0, pi/2, num_points);  % 角度范围 (1/4 圆柱)
z_cyl = linspace(z_min, z_max, num_points)';  % 高度方向划分

[Theta, Z_cyl] = meshgrid(theta, z_cyl);
R = radius * ones(size(Theta));  % 设定固定半径

% **转换为直角坐标**
X_cyl = R .* cos(Theta);
Y_cyl = R .* sin(Theta);

hold on;

% **绘制 1/4 透明圆柱**
surf(X_cyl, Y_cyl, Z_cyl, 'FaceColor', cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value);

% **封闭底部（1/4圆形，透明）**
theta_bottom = linspace(0, pi/2, num_points);  % 1/4 圆的角度
X_bottom = radius * cos(theta_bottom);
Y_bottom = radius * sin(theta_bottom);
Z_bottom = z_min * ones(size(X_bottom));  % 设置底部Z值
fill3([0 X_bottom], [0 Y_bottom], [z_min Z_bottom], cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value); % 底部透明封闭

% **封闭顶部（1/4圆形，透明）**
Z_top = z_max * ones(size(X_bottom));  % 设置顶部Z值
fill3([0 X_bottom], [0 Y_bottom], [z_max Z_top], cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value); % 顶部透明封闭

% **在圆柱顶面中心添加 "rod" 文字**
text(radius/2, radius/2, z_max, 'Fuel rod', 'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'r', 'HorizontalAlignment', 'center');

hold off;







% **Compute the area of the juxing_area + P2_area
dx = diff(X_combined, 1, 2); % Spacing in the X direction
dy = diff(Y_combined, 1, 1); % Spacing in the Y direction

% **Ensure dx and dy do not contain NaN values**
dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

% **Calculate the area of juxing_area + P2_area
cell_area_combined = dx(2:end, :) .* dy(:, 2:end); 

% **Calculate the area of the triangular region**
tri_area = 0.5 * W * W;  % Directly use the formula A = 1/2 * W^2

% **Compute valid flow field values in the juxing_area + P2_area
valid_mask_combined = ~isnan(F_combined(2:end, 2:end));
F_valid_combined = F_combined(2:end, 2:end);
F_valid_combined(~valid_mask_combined) = 0;
cell_area_combined(~valid_mask_combined) = 0;

% **Calculate the actual flow field area of the juxing_area + P2_area
actual_area_combined = sum(cell_area_combined(:));

% **Compute the average flow field value in the triangular region**
F_tri_avg = mean(F_tri(:), 'omitnan'); % Directly use the mean value of the triangular region

% **Total interpolated flow field area**
actual_area_total = actual_area_combined + tri_area;

% **Compute the area-weighted average flow field value**
if actual_area_total == 0
    area_weighted_avg = NaN;
else
    area_weighted_avg = (sum(F_valid_combined(:) .* cell_area_combined(:)) + F_tri_avg * tri_area) / actual_area_total;
end

% **Compute the Paranjape value**
Paranjape = 0.68 * f_0 + 0.22 * f_g;

% **Output results**
fprintf('juxing_area + P2_area (after RBF): %.4f\n', actual_area_combined);
fprintf('Triangular region area: %.4f\n', tri_area);
fprintf('juxing_area + P2_area + Triangular_area (after RBF): %.4f\n', actual_area_total);
fprintf('Area-weighted average flow field value: %.4f\n', area_weighted_avg);
fprintf('Paranjape area-weighted average flow field value: %.4f\n', Paranjape);


% **RBF 训练效果评估**
% **拆分数据集**
num_samples = size(X_train, 1);
rand_indices = randperm(num_samples);  % 随机打乱索引
train_size = round(0.8 * num_samples); % 80% 训练，20% 测试

train_indices = rand_indices(1:train_size);
test_indices = rand_indices(train_size+1:end);

X_train_sub = X_train(train_indices, :);
F_train_sub = F_train(train_indices);

X_test = X_train(test_indices, :);
F_test = F_train(test_indices);

% **重新计算 RBF 权重 (仅使用 80% 训练数据)**
D_train = pdist2(X_train_sub, X_train_sub);
sigma = 0.15 * median(D_train(:));  % 计算 sigma
rbf_func = @(r) exp(-(r / sigma).^2);
Phi_train = rbf_func(D_train);
lambda = 1e-6;  
Phi_train = Phi_train + lambda * eye(size(Phi_train));
weights = Phi_train \ F_train_sub;

% **对测试数据集进行预测**
D_test = pdist2(X_test, X_train_sub);
Phi_test = rbf_func(D_test);
F_pred = Phi_test * weights; % 预测流场值

% **计算误差**
mse = mean((F_test - F_pred).^2);  % 均方误差 (MSE)
rmse = sqrt(mse);  % 均方根误差 (RMSE)
ss_total = sum((F_test - mean(F_test)).^2);  % 总体平方和
ss_residual = sum((F_test - F_pred).^2);  % 残差平方和
r2 = 1 - (ss_residual / ss_total);  % R²（决定系数）

% **输出误差指标**
fprintf('RBF 训练效果评估:\n');
fprintf('均方误差 (MSE): %.6f\n', mse);
fprintf('均方根误差 (RMSE): %.6f\n', rmse);
fprintf('决定系数 (R²): %.4f\n', r2);

% **绘制预测值与真实值的对比**
% 归一化 F_test 和 F_pred（Min-Max 归一化）
F_test_norm = (F_test - min(F_test)) / (max(F_test) - min(F_test));
F_pred_norm = (F_pred - min(F_test)) / (max(F_test) - min(F_test));

figure(8);
scatter(F_test_norm, F_pred_norm, 30, 'b', 'filled'); hold on;
plot([min(F_test_norm), max(F_test_norm)], [min(F_test_norm), max(F_test_norm)], 'r--', 'LineWidth', 1.5);
title('80% for training, 20% for predicting');

% 设置轴标签
xlabel('\it F\rm_{test}, -', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('\it F\rm_{pred}, -', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');

% 设置 legend 并调整符号大小
h = legend('Predicted values', 'True values', 'Location', 'northwest');
set(h, 'Box', 'off');                
set(h, 'FontName', 'Times New Roman'); 
set(h, 'FontSize', 14); 
set(gcf, 'Renderer', 'painters'); % 使用矢量渲染
h.ItemTokenSize = [20, 10]; % **放大 legend 符号**

daspect([1 1 0.5]); % X and Y axes equal, Z scaled dynamically

% **设置坐标轴格式，确保所有刻度可见**
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
set(gca, 'LineWidth', 2);
box on; % **确保所有框线显示**

% **确保所有刻度都显示**
xticks(min(F_test_norm):0.2:max(F_test_norm)); % X 轴刻度
yticks(min(F_test_norm):0.2:max(F_test_norm)); % Y 轴刻度

% **设置轴颜色确保所有刻度都可见**
set(gca, 'XColor', 'k', 'YColor', 'k'); % 确保 X、Y 轴的刻度线显示
set(gca, 'TickDir', 'in');  % 让刻度向内

% 设置 x 和 y 轴范围
xlim([0 1]);
ylim([0 1]);

% 找到图的中心
x_center = (min(F_test_norm) + max(F_test_norm)) / 2;
y_center = (min(F_pred_norm) + max(F_pred_norm)) / 2;

% 计算新的 R²、MSE 和 RMSE
r2_norm = 1 - sum((F_test_norm - F_pred_norm).^2) / sum((F_test_norm - mean(F_test_norm)).^2);
mse_norm = mean((F_test_norm - F_pred_norm).^2);
rmse_norm = sqrt(mse_norm);

% 在图中间标注 R² 值
text(x_center-0.2, y_center+0.1, sprintf('R^{2} = %.3f\nMSE = %.3f\nRMSE = %.3f', r2_norm, mse_norm, rmse_norm), ...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'center', 'BackgroundColor', 'w','color','r');

hold off;
fprintf('Min/Max of F_test (normalized): %.4f / %.4f\n', min(F_test_norm), max(F_test_norm));
fprintf('Min/Max of F_pred (normalized): %.4f / %.4f\n', min(F_pred_norm), max(F_pred_norm));

