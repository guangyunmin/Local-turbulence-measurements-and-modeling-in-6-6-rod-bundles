clc;
clear;
close all;

exp_data = importdata('2-18.txt');
exp_data = exp_data.data;


% Parameter settings
Re_t_values = [180, 395, 590, 1000, 1.35*0.4*1*10^6/log(1*10^6)/2];  % Different turbulent Reynolds numbers
colors = {'r', 'k', 'b', 'g', 'm'};  % Define colors as cell array

% Iterate over Reynolds numbers
for i = 1:length(Re_t_values)
    Re_t = Re_t_values(i);
    k = @(y) 0.4 - (0.4 - 0.14) .* (y);
    
    % Define L(y)
    L_1 = @(y) (0.14 - 0.08 .* (y).^2 - 0.06 * (y).^4);
    L_2 = @(y) L_1(y) .* (1 - exp((y - 1) .* Re_t / 26));
    L_3 = @(y) L_1(y) .* (1 - exp((y - 1) .* Re_t / 26)) ...
        ./ sqrt(max(1 - exp((y-1) .* Re_t .* 0.26), 1e-10)); 
    L_4 = @(y) 0.4 .* (1 - y) .* (1 - exp((y - 1) .* Re_t / 26));
    L_5 = @(y) (0.4 .* (1 - y));

    % Define differential equations
    odefun = @(L) @(y, u) (-2 * y * Re_t) ./ (1 + sqrt(1 + 4 * (L(y)).^2 * Re_t^2 .* y));

    % Solve ODEs
    yspan = 1:-0.001:0;  % Range of y
    u0 = 0;              % Initial condition
    [y_1, u_1] = ode45(odefun(L_1), yspan, u0);
    [y_2, u_2] = ode45(odefun(L_2), yspan, u0);
    [y_3, u_3] = ode45(odefun(L_3), yspan, u0);
    [y_4, u_4] = ode45(odefun(L_4), yspan, u0);
    [y_5, u_5] = ode45(odefun(L_5), yspan, u0);

    % Plot results
    figure(i);
    hold on;
    plot(y_1, u_1/max(u_1), 'r', 'DisplayName', 'Van Driest model', 'LineWidth', 2);
    plot(y_2, u_2/max(u_2), 'k', 'DisplayName', 'Cebeci and Bradshaw model', 'LineWidth', 2);
    plot(y_3, u_3/max(u_3), 'b', 'DisplayName', 'Hanna model', 'LineWidth', 2);
    plot(y_4, u_4/max(u_4), 'g', 'DisplayName', 'Nikuradse model', 'LineWidth', 2);
    plot(y_5, u_5/max(u_5), 'm', 'DisplayName', 'Prandtl model', 'LineWidth', 2);

    % DNS data
    plot(exp_data(:,1), exp_data(:,2)/max(exp_data(:,2)), 'y', 'DisplayName', 'Exp data', 'LineWidth', 2);


    xlabel('r/R');
    ylabel('Velocity');
    title(['Velocity Profile Fit (Re_t = ' num2str(Re_t) ')']);
    legend('show', 'box', 'off', 'location', 'best');
    grid on;
    hold off;
end
