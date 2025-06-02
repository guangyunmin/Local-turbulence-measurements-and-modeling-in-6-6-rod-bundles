clc; clear; close all;
load('Data_central.mat')
% **拟合 Fina(F0, FG) 的显式函数表达式**
% 选择一个合适的拟合方法，确保函数形式易解析
% `poly22` 代表 F0 为二次项, FG 为三次项的多项式拟合

% 将 FO 和 FG 展平，并存储 Fina_avg 中的值
FO_flat = FO(:);
FG_flat = FG(:);
Fina_flat = Fina_avg(:);
Fina_fit = fit([FO_flat, FG_flat], Fina_flat, 'poly22'); 

% **获取拟合函数的解析表达式**
Fina_expr = formula(Fina_fit);

% **获取系数**
coeffs = coeffvalues(Fina_fit);
Fina_coeffs = coeffnames(Fina_fit);

% **打印显式函数关系式**
fprintf('Fina(F0, FG) = %s\n', Fina_expr);

% **替换表达式中的系数变量**
for i = 1:length(Fina_coeffs)
    Fina_expr = strrep(Fina_expr, Fina_coeffs{i}, sprintf('%.6f', coeffs(i)));
end

fprintf('具体拟合关系式:\nFina(F0, FG) = %s\n', Fina_expr);


% **测试输出**
test_F0 = 1;
test_FG = 0.97;
Fina_test_value = Fina_fit(test_F0, test_FG);
fprintf('Fina(%.2f, %.2f) = %.6f\n', test_F0, test_FG, Fina_test_value);


% **绘制拟合曲面和原始数据点在同一图**
figure(9);
hold on;

% **绘制拟合曲面**
[X_fit, Y_fit] = meshgrid(linspace(min(FO(:)), max(FO(:)), 50), ...
                          linspace(min(FG(:)), max(FG(:)), 50));
Z_fit = Fina_fit(X_fit, Y_fit); % 计算拟合值

surf(X_fit, Y_fit, Z_fit, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); % 透明度 0.8
colormap('jet');
colorbar;

% **绘制原始数据散点**
scatter3(FO_flat, FG_flat, Fina_flat, 50, 'r', 'filled'); % 红色原始数据点

% **设置坐标轴**
xlabel('F0', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('FG', 'FontSize', 14, 'FontName', 'Times New Roman');
zlabel('Fina(F0, FG)', 'FontSize', 14, 'FontName', 'Times New Roman');
h_legend=legend('Fitted function', 'Data points', 'Location', 'best','box','off');
set(h_legend, 'Position', [0.4, 0.7, 0.1, 0.1])

% **Set labels with Times New Roman Italic**
xlabel('\it X\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('\it Y\rm, mm', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
zlabel('\it Z\rm, m/s', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter', 'tex');
title('1 / 8 center subchannel', 'FontName', 'Times New Roman', 'FontSize', 14);

% **Ensure X and Y axes have equal scaling but leave Z free**
% daspect([1 1 1]); % X and Y axes equal, Z scaled dynamically

% **Set tick intervals to 1**
% xticks(min(FO(:)):1:max(FO(:)));
% yticks(min(FG(:)):1:max(FG(:)));
% % % zticks(min( F_final):1:max( F_final));
% xlim([0 max(FO(:))]);
% ylim([0 max(FG(:))]);
% zlim([0 max(F_final)]);
% **Set grid and view**
grid on;
view(30, 40);
% **设置 X, Y, Z 轴的线条更粗**
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');
set(gcf, 'Renderer', 'painters'); % 使用矢量渲染
hold off;

X_fit=X_fit(:);
Y_fit=Y_fit(:);
Z_fit=Z_fit(:);
XX=[FO_flat,FG_flat,Fina_flat];
XX_fit=[X_fit,Y_fit,Z_fit];

