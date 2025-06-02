clc; 
clear; 
close all;

% =========================
% 🎯 Calculate Subchannel Area
% =========================
P = 16.7;
D = 10;
W = 1.665; 
Across = (100^2-(36*pi*(D/2)^2))/4;


% **读取并处理流场数据**
% 你的数据点
data = [
    0	0	0.533946519366532;
    0	8.3333	0.540123880831321;
    0	25	0.497679098793165;
    0	16.6667	0.523677633289183;
    -8.3333	16.6667	0.509968965886419;
    -16.6667	16.6667	0.506793718979753;
    -16.6667	25	0.500335630242869;
    8.3333	16.6667	0.484868880172278;
    16.6667	16.6667	0.484061807625539;
    16.6667	25	0.483380257124396;
    16.6667	41.6667	0.459234793616479;
    0	33.3333	0.516032847653823;
    -8.3333	33.3333	0.477232057773877;
    -16.6667	33.3333	0.506366631523833;
    -25	33.3333	0.500174560697931;
    -33.3333	33.3333	0.505408119825444;
    8.3333	33.3333	0.507409904393652;
    16.6667	33.3333	0.513236253517183;
    25	33.3333	0.499815048415271;
    33.3333	33.3333	0.503072120375269;
    33.3333	41.6667	0.473501042710423;
    0	41.6667	0.504342099079284;
    -16.6667	41.6667	0.475350621303607;
    -33.3333	41.6667	0.489053452735998;
    0	48.335	0.502522658192896;
    -8.3333	48.335	0.460194055335227;
    -16.6667	48.335	0.504935714639818;
    -25	48.335	0.487702315924416;
    -33.3333	48.335	0.503226004597474;
    -41.6667	48.335	0.447279853218646;
    8.3333	48.335	0.463032397411362;
    16.6667	48.335	0.483545482656568;
    25	48.335	0.478491259905534;
    33.3333	48.335	0.482711710767044;
    41.6667	48.335	0.438946259441051;
    -48.335	48.335	0.457684461816915;
    48.335	48.335	0.462299016993349
];
save('data1.mat','data')
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


figure;
scatter(unique_points(:,1), unique_points(:,2), 50, v_avg, 'filled');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('完整的流场数据分布（去重后）');
grid on;
axis equal;
colormap('jet')




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




