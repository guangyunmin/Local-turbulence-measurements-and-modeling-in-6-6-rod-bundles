clc; clear; close all;

% 定义参数
Re_tau = 180;  % 例子值，可修改


y30 = (1 - sqrt(1 + 2 * Re_tau^2)) / Re_tau - sqrt(2)/(4*Re_tau^2);
% y30 =-sqrt(2) + 1/Re_tau;
(49*Re_tau^3*(exp(-Re_tau/26) - 1)^2)/1250;

% 定义求解区间（修正 r 递减问题）
r_span = 0:0.005:1;
y0 = [33; 0; y30]; % 初始条件 [u1, u2, u3, u4]

% 调用 ODE45 求解
[r_values, Y] = ode45(@(r, y) odefunc(r, y, Re_tau), r_span, y0);

% 画图
figure;
plot(r_values, Y(:, 1), 'b', 'LineWidth', 2);
xlabel('$\tilde{r}$', 'Interpreter', 'latex');
ylabel('$\tilde{u}$', 'Interpreter', 'latex');
title('Solution using ODE45', 'Interpreter', 'latex');
grid on;
save('my_model.mat','r_values','Y')
% ------------------- 函数定义 ------------------- %
% 定义 ODE 系统
function dydr = odefunc(r, y, Re_tau)
    y1 = y(1); 
    y2 = y(2); 
    y3 = y(3); 

    % 计算 lambda
    lambda = (0.14 - 0.08 .* r.^2 - 0.06 .* r.^4) .* (1 - exp((r - 1) .* Re_tau / 26));

    % 计算 alpha, beta, chi
    alpha = Re_tau * lambda^2;
    beta = (1/2) * Re_tau * lambda^3;
    chi = (1/6) * Re_tau * lambda^4;

    % 避免 chi 过小导致除零问题
    if abs(chi) < 1e-18
        chi = 1e-17;
    end

    % 避免 y2 过小导致数值不稳定
    if abs(y2) < 1e-18
        y2 = 1e-17;
    end

    % 计算 dy3/dr
    dydr3 = (1 - alpha * y2 - beta * y3 + Re_tau * r / y2) / chi;

    % 返回微分方程组
    dydr = [y2; y3; dydr3];
end
