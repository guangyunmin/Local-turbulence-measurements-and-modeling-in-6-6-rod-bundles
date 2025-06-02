clc; clear; close all;

%% 参数定义
R = 1;
Re_tau = 180;
gamma = 26;
kappa = 0.4;

% 网格
Nx = 200; Ny = 200;
x = linspace(-R, R, Nx);
y = linspace(-R, R, Ny);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);
r_tilde = 1 - r / R;
mask = r <= R;

% 通用代数项
base = (1 - 4/7 * (1 - r_tilde).^2 - 3/7 * (1 - r_tilde).^4);

%% 模型 A：基础代数模型
lambdaA = 0.14 * base;
lambdaA(~mask) = NaN;

%% 模型 B：代数 × 指数抑制
lambdaB = lambdaA .* (1 - exp(-Re_tau * r_tilde / gamma));
lambdaB(~mask) = NaN;

%% 模型 C：Wall函数归一化指数修正
denominator = sqrt(1 - exp(-0.26 * Re_tau* r_tilde));
lambdaC = lambdaA .* (1 - exp(-Re_tau * r_tilde / gamma)) ./ denominator;
lambdaC(~mask) = NaN;

%% 模型 D：线性 von Kármán 模型
lambdaD = kappa * r_tilde;
lambdaD(~mask) = NaN;

%% 模型 E：指数型模型（不含代数项）
lambdaE = 0.4 * r_tilde .* (1 - exp(-Re_tau * r_tilde / gamma));
lambdaE(~mask) = NaN;

%% 可视化：2x3 subplot
figure;
models = {lambdaA, lambdaB, lambdaC, lambdaD, lambdaE};
titles = {'A: 基础代数模型', ...
          'B: 指数抑制代数模型', ...
          'C: 归一化 Wall函数模型', ...
          'D: 线性 Kármán 模型', ...
          'E: 纯指数模型'};

for i = 1:5
    figure(i);
    surf(X, Y, models{i}, 'EdgeColor', 'none');
    colormap(jet);
%     zlim([0 0.14]);
    title(titles{i}, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('X', 'FontName', 'Times New Roman');
    ylabel('Y', 'FontName', 'Times New Roman');
    zlabel('\lambda', 'FontName', 'Times New Roman');
    view(45, 30);
    colorbar;
    
end
XX=X(:);
YY=Y(:);
ZZ=models{i}(:);
% 样式统一
set(findall(gcf,'Type','axes'), ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.5);
