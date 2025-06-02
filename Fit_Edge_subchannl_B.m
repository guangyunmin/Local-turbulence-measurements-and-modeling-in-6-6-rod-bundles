clc; clear; close all;

Fina_avg=zeros(size( 0.4:0.1:1.2,2),size( 0.5:0.1:1.2,2));
% f_g1 = [0.46/0.5, 0.64/0.70,0.74/0.82,0.88/0.98,1.12/1.20]
f_g1 = 0.4:0.1:1.2;
f_01 = 0.5:0.1:1.2;
[FO,FG] = meshgrid(f_01,f_g1);
for kk=1:size(f_g1,2)
    for kj=1:size(f_01,2)
        f_g = f_g1(1,kk);
        f_0 = f_01(1,kj);
        
        % **Parameter Definition**
        P = 16.7;   % Structure size parameter
        D = 10;     % Diameter of the circular region
        m = 7;      % Power exponent
        W = 1.665;
        R0 = D / 2;   % Initial radius
        Rnew = 0.72*R0;
        
        % **Create polar coordinate grid**
        theta = linspace(0, pi / 4, 50);
        r = linspace(0, (P/2-W) * sqrt(2), 50);
        [THETA, R] = meshgrid(theta, r);
        
        % **Compute flow field**（图2需要条件这
        sec_theta = 1 ./ cos(THETA);
        f_c = f_g + (f_0 - f_g) .* (1 - ((sqrt(2) - sec_theta) / (sqrt(2) - 1)).^m);
        F = f_c .* (1 - (1 - (R - Rnew ) ./ ((P / 2-W) .* sec_theta -Rnew )).^m);
        F(R <= Rnew ) = 0;
        
        % **Convert to Cartesian coordinates**
        [X, Y] = pol2cart(THETA, R);
        X_data = [X(:), Y(:)];
        F_data = F(:);
        
        % **Compute P2 structured grid**
        theta_hou = atan((P/2-W) / (P/2));
        arc_theta_hou = linspace(theta_hou, 0, 30);
        arc_theta_qian = linspace(pi/4, 0, 30);
        R_qian = 0.72*R0;
        R_hou = R0;
        arc_x_qian = R_qian * cos(arc_theta_qian);
        arc_y_qian = R_qian * sin(arc_theta_qian);
        arc_x_hou = R_hou * cos(arc_theta_hou)-W;
        arc_y_hou = R_hou * sin(arc_theta_hou);
        arc_points_hou = [arc_x_hou', arc_y_hou'];
        arc_points_qian = [arc_x_qian', arc_y_qian'];
        B = [P/2-W, 0];
        C = [P/2-W, P/2-W];
        D_hou = [R_hou * cos(theta_hou)-W , R_hou * sin(theta_hou)];
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
        juxing_area = 0 * (P / 2);
        
        % **Output the areas of P1 and P2**
        fprintf('P1 structured grid area: %.4f\n', P1_area);
        fprintf('P2 structured grid area: %.4f\n', P2_area);
        fprintf('Rectangle area: %.4f\n', juxing_area);
        fprintf('Total area before interpolation: %.4f\n', juxing_area + P2_area);
        
        
        % **Create P2 structured grid**
        P2_x = linspace(min(P2_points(:,1)), max(P2_points(:,1)), 300);
        P2_y = linspace(min(P2_points(:,2)), max(P2_points(:,2)), 300);
        [P2_X, P2_Y] = meshgrid(P2_x, P2_y);
        
        [in, ~] = inpolygon(P2_X, P2_Y, P2_points(:,1), P2_points(:,2));
        P2_X(~in) = NaN;
        P2_Y(~in) = NaN;
        
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
        sigma = 0.15*median(D_train(:));
        rbf_func = @(r) exp(-(r / sigma).^2);
        Phi_train = rbf_func(D_train);
        lambda = 1e-6;
        Phi_train = Phi_train + lambda * eye(size(Phi_train));
        weights = Phi_train \ F_train;
        
        D_target = pdist2(X_target, X_train);
        Phi_target = rbf_func(D_target);
        
        F_P2_grid = NaN(size(P2_X));
        F_P2_grid(valid_target_idx) = Phi_target * weights;
        
        
        
        % **Merge P2 grid and right rectangular grid**
        X_combined = [P2_X];
        Y_combined = [P2_Y];
        F_combined = [F_P2_grid];
        
        F_combined(isnan(F_combined)) = 0; % Avoid NaN affecting visualization
        
        % **Figure 2: Flow field before interpolation (Cartesian coordinates)**
        figure(2);
        surf(X, Y, F, 'EdgeColor', 'none');
        colormap('jet');
        colorbar;
        xlabel('X'); ylabel('Y'); zlabel('F');
        title('Flow field before interpolation (Cartesian coordinates)');
        
        % **Figure 3: Interpolated complete flow field**
        figure(3);
        surf(P2_X, P2_Y, F_P2_grid, 'EdgeColor', 'none');
        colormap('jet');
        colorbar;
        xlabel('X'); ylabel('Y'); zlabel('F');
        title('Interpolated flow field');
        
        % **Figure 4: Grid shape comparison before and after interpolation (Cartesian coordinates)**
        figure(4);
        hold on;
        scatter(X_data(:,1), X_data(:,2), 10, 'b', 'filled');
        scatter(P2_X(:), P2_Y(:), 10, 'r', 'filled');
        
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
        
        
        
        % **Figure 5: Merged flow field**
        
        
        figure(5);
        surf(X_combined, Y_combined, F_combined, 'EdgeColor', 'none');
        colormap(jet);
        colorbar;
        
        % **Set labels with Times New Roman Italic**
        xlabel('\it X\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
        ylabel('\it Y\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
        zlabel('\it Z\rm, -', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
        title('1 / 4 edge subchannel (type B)', 'FontName', 'Times New Roman', 'FontSize', 14);
        
        % **Ensure X and Y axes have equal scaling but leave Z free**
        daspect([1 1 0.5]); % X and Y axes equal, Z scaled dynamically
        
        % **Set tick intervals to 1**
        xticks(0:1:max(X_combined(:)));
        yticks(min(Y_combined(:)):1:max(Y_combined(:)));
        % % zticks(min( F_final):1:max( F_final));
        xlim([-W 7]);
        ylim([0 max(Y_combined(:))]);
        % zlim([0 max(F_final)]);
        % **Set grid and view**
        grid on;
        view(30, 40);
        % **设置 X, Y, Z 轴的线条更粗**
        set(gca, 'LineWidth', 2);
        set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
        
        
        
        
        
        % **Find the maximum value of F_final and its index**
        [max_F, max_idx] = max(F_combined(:)); % Find the maximum value
        max_X = X_combined(max_idx); % Get the corresponding X coordinate
        max_Y = Y_combined(max_idx); % Get the corresponding Y coordinate
        hold on;
        plot3(max_X, max_Y, max_F, 'bo', 'MarkerSize', 14, 'MarkerFaceColor', 'b'); % Mark maximum value
        text(max_X-W/2, max_Y, max_F+0.4,  sprintf('\\it f\\rm_{c}: %.2f', max_F), ...
            'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'r', 'Interpreter', 'tex');
        scatter3(D/2+W, 0, 0.98, 200, [0.5, 1, 1], 'filled');
        text(D/2+W, -2.4, 0.98+0.4,sprintf('\\it f\\rm_{g}: %.2f', 0.98), ...
            'FontSize', 14, 'FontName', 'Times New Roman', ...
            'Interpreter', 'tex', 'Color', 'r');
        
        % **Add a cylinder (fuel rod) centered at (0,0)**
        
        
        
        % **定义圆柱参数**
        x_center = -W;  % 圆柱中心 x 变为 -W
        y_center = 0;   % 圆柱中心 y 不变
        z_min = min(F_combined(:)); % 圆柱底部
        z_max = max(F_combined(:)); % 圆柱顶部
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
        X_cyl = R .* cos(Theta) - W; % 调整中心到 (-W, 0)
        Y_cyl = R .* sin(Theta);
        
        hold on;
        
        % **绘制 1/4 透明圆柱**
        surf(X_cyl, Y_cyl, Z_cyl, 'FaceColor', cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value);
        
        % **封闭底部（1/4圆形，透明）**
        theta_bottom = linspace(0, pi/2, num_points);  % 1/4 圆的角度
        X_bottom = radius * cos(theta_bottom) - W; % 调整 X 坐标
        Y_bottom = radius * sin(theta_bottom);
        Z_bottom = z_min * ones(size(X_bottom));  % 设置底部Z值
        fill3([x_center X_bottom], [y_center Y_bottom], [z_min Z_bottom], cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value); % 底部透明封闭
        
        % **封闭顶部（1/4圆形，透明）**
        Z_top = z_max * ones(size(X_bottom));  % 设置顶部Z值
        fill3([x_center X_bottom], [y_center Y_bottom], [z_max Z_top], cyl_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value); % 顶部透明封闭
        
        % **在圆柱顶面中心添加 "rod" 文字**
        text(x_center + radius/2, y_center + radius/2, z_max+ 0.2, 'Fuel rod', ...
            'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'r', ...
            'HorizontalAlignment', 'center');
        
        
        % set(gcf, 'Renderer', 'painters'); % 使用矢量渲染
        % print('-clipboard', '-dmeta');     % 复制到剪贴板为矢量格式
        hold off;
        
        
        
        
        
        
        
        
        
        
        
        % **Compute grid spacing in X and Y directions**
        dx = diff(X_combined, 1, 2); % X (size: rows, cols-1)
        dy = diff(Y_combined, 1, 1); % Y (size: rows-1, cols)
        
        % **Ensure dx and dy do not contain NaN values**
        dx(isnan(dx)) = 0;
        dy(isnan(dy)) = 0;
        
        % **Compute the area of grid cells**
        dx_avg = 0.5 * (dx(1:end-1, :) + dx(2:end, :));
        dy_avg = 0.5 * (dy(:, 1:end-1) + dy(:, 2:end));
        
        % **Ensure area calculation does not contain NaN values**
        cell_area = dx_avg .* dy_avg;
        cell_area(isnan(cell_area)) = 0; % Avoid NaN affecting the calculation
        
        % **Debugging: Print the minimum and maximum grid cell area**
        fprintf('Debug Info: Minimum computed grid cell area: %.4f, Maximum computed grid cell area: %.4f\n', min(cell_area(:)), max(cell_area(:)));
        
        % **Filter valid regions of F_combined (remove NaN values)**
        valid_mask = ~isnan(F_combined(2:end, 2:end));
        
        % **Ensure NaN regions do not affect the calculation**
        cell_area(~valid_mask) = 0;
        F_valid = F_combined(2:end, 2:end);
        F_valid(~valid_mask) = 0;
        
        % **Compute the actual 2D flow field area**
        actual_area = sum(cell_area(:));
        
        % **Prevent division by zero**
        if actual_area == 0
            area_weighted_avg = NaN;
        else
            area_weighted_avg = sum(F_valid(:) .* cell_area(:)) / actual_area;
        end
        
        % **Final output**
        Paranjape=0.68*f_0+0.22*f_g;
        fprintf('Interpolated flow field area (2D calculation): %.3f\n', actual_area);
        fprintf('Flow field area weighted average: %.3f\n', area_weighted_avg);
        fprintf('Paranjape area weighted average flow field value: %.3f\n', Paranjape);
        
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
        
        Fina_avg(kk,kj)=area_weighted_avg ;
    end
end

save('Data_edge_B.mat','FO','FG','Fina_avg')