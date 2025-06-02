clc; clear; close all;
% ================== 计算混合长度 (lm) 的分布 ==================
% Set parameters
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
decay_factor = Re_t_values(3) / C4;
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

 
% 通道尺寸和幂律参数
W = 6 * P;  % 这里选择流道宽度为 6P，与混合长度网格对齐
H = 6 * P;  % 同理

% 为了配合混合长度计算的 X, Y 网格，重新归一化坐标
r_x = abs(X) / (W/2);
r_y = abs(Y) / (H/2);

% 通道幂律分布（使用均匀最大值归一化）
% U_channel_x = 0.14*(1-4/7*(r_x).^2-3/7*(r_x).^4).*(1-exp(-Re_t_values(1).*(1-r_x)/C4))/(W/2);
% U_channel_y = 0.14*(1-4/7*(r_y).^2-3/7*(r_y).^4).*(1-exp(-Re_t_values(1).*(1-r_y)/C4))/(W/2);
% 
% U_channel =( U_channel_x.^(-20) + U_channel_y.^(-20)).^(-1/20);


U_channel_x = (1-4/7*(r_x).^2-3/7*(r_x).^4).*(1-exp(-Re_t_values(1).*(1-r_x)/C4));
U_channel_y = (1-4/7*(r_y).^2-3/7*(r_y).^4).*(1-exp(-Re_t_values(1).*(1-r_y)/C4));

U_channel =0.14*( U_channel_x.* U_channel_y);

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


% 可视化混合长度 lm 分布
figure;
surf(X, Y, lm.* U_scaled,  'EdgeColor', 'none'); 
AAA=lm.* U_scaled;
XX=X(:);
YY=Y(:);
ZZ=(AAA(:));
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


% 可视化混合长度 lm 分布
figure;
surf(X, Y, U_channel,  'EdgeColor', 'none'); 

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
XY_target=unique_points;

% 提取 X 与 Y 网格中的间距（确保插值有效）
x_range = linspace(-3*P, 3*P, Nx)*1000;
y_range = linspace(-3*P, 3*P, Ny)*1000;

% 插值计算目标点处的 AAA 值
AAA_values = interp2(X*1000, Y*1000, AAA, x, y, 'linear');


U_s=[0	0	0.0114441253771623
0	8.3333	0.0103140469272618
0	25	0.0120593387656858
		
		
		
		
0	16.6666666666667	0.0108340250978894
-8.33333333333333	16.6666666666667	0.0110457086277028
-16.6666666666667	16.6666666666667	0.0112781601682003
-16.6667	25	0.0112893965673789
		
8.33333333333333	16.6666666666667	0.0114867355060166
16.6666666666667	16.6666666666667	0.0113472687852212
16.6667	25	0.0114407943418863
16.6667	41.6667	0.0117987658994313
0	33.3333333333333	0.0115213830309157
-8.33333333333333	33.3333333333333	0.0116821077822736
-16.6666666666667	33.3333333333333	0.0104729025885856
-25	33.3333333333333	0.0119847628928235
-33.3333333333333	33.3333333333333	0.0122860002963452
		
8.33333333333333	33.3333333333333	0.010899895069688
16.6666666666667	33.3333333333333	0.0102661702031115
25	33.3333333333333	0.0105453935654542
33.3333333333333	33.3333333333333	0.0106195924084627
33.3333	41.6667	0.0109598157369238
		
		
		
		
		
		
		
		
		
		
		
		
		
		
0	41.6667	0.0119847631972972
-16.6667	41.6667	0.01181
-33.3333	41.6667	0.0105366022377694
0	48.335	0.0121285218874169
-8.33333333333333	48.335	0.0130123917201042
-16.6666666666667	48.335	0.0129309134644813
-25	48.335	0.0129701675881712
-33.3333333333333	48.335	0.0113690215930002
-41.6666666666667	48.335	0.0129568524818746
		
8.33333333333333	48.335	0.012312718480826
16.6666666666667	48.335	0.0117337485442269
25	48.335	0.0118821639291898
33.3333333333333	48.335	0.0118204346553687
41.6666666666667	48.335	0.0129090184807932
-48.335	48.335	0.0126667224416452
48.335	48.335	0.0120656730347284];
% 载入数据

Eddy_frequcy = U_s(:,3)./ (AAA_values * 0.1);
min(Eddy_frequcy(:))
max(Eddy_frequcy(:))

