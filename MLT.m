clc; 
clear; 
close all;

% Parameter settings
Re_t_values = [180, 395, 590, 1000];  % Different turbulent Reynolds numbers
colors = ['b', 'k', 'r'];  % Add black for DNS data
u_DNS_180 = importdata("Re_180.txt");
u_DNS_180 = u_DNS_180.data;
u_DNS_395 = importdata("Re_395.txt");
u_DNS_395 = u_DNS_395.data;
u_DNS_590 = importdata("Re_590.txt");
u_DNS_590 = u_DNS_590.data;
u_DNS_1000 = importdata("Re_1000.txt");
u_DNS_1000 = u_DNS_1000.data;
u_DNS = [u_DNS_180, u_DNS_395, u_DNS_590, u_DNS_1000];
save('u_DNS.mat','u_DNS')

for i = 1:length(Re_t_values)
    Re_t = Re_t_values(i);
    
    % Define L(y)
    L_1 = @(y) (0.14 - 0.08 .* (y).^2 - 0.06 * (y).^4);
    L_2 = @(y) (0.14 - 0.08 .* (y).^2 - 0.06 * (y).^4).*(1 - exp((y-1)*Re_t/26));

    % Define differential equations
    odefun_1 = @(y, u) (-2 * y * Re_t) ./ (1 + sqrt(1 + 4 * (L_1(y)).^2 * Re_t^2 .* y));
    odefun_2 = @(y, u) (-2 * y * Re_t) ./ (1 + sqrt(1 + 4 * (L_2(y)).^2 * Re_t^2 .* y));
    
    % Set the solution range
    yspan = [1 0];  % Set the range of y
    u0 = 0;         % Set initial conditions
    
    % Call ode45 to solve
    [y_1, u_1] = ode45(odefun_1, yspan, u0);
    [y_2, u_2] = ode45(odefun_2, yspan, u0);
    
    % Normalize u
    u_norm_1 = u_1 / max(u_1);
    u_norm_2 = u_2 / max(u_2);
    
    % Plot an independent figure for each Re_t
    figure(i);
    hold on;
    plot(u_norm_1 , y_1, [colors(1) '-'], 'LineWidth', 2, 'DisplayName', 'No Damping');
    plot(u_norm_2 , y_2, [colors(2) '--'], 'LineWidth', 2, 'DisplayName', 'Damping');
    plot(u_DNS(:, 2*i-1) , u_DNS(:, 2*i), [colors(3) 'o'], 'LineWidth', 2, 'DisplayName', 'DNS');
    hold off;
    
    ylabel('$\mathit{y}/\mathit{R}$', 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('$\bar{u} / \bar{u_c}$', 'Interpreter', 'latex', 'FontSize', 16);
    title(['Turbulent Mean Velocity Profile (Re_t = ' num2str(Re_t) ')'], 'FontSize', 12);
    xlim([0,1]);
    legend('show', 'Location', 'best');
    grid on;
end
