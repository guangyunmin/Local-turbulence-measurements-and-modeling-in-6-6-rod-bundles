% 定义 R_hat
R_hat = 10;

% 计算原始公式的值
result_original = (1 - sqrt(1 + 2 * R_hat^2)) / R_hat;

% 计算渐近展开公式的值
approx_result = -sqrt(2) + 1/R_hat;

% 显示结果
fprintf('原始公式计算结果: %.10f\n', result_original);
fprintf('渐近展开计算结果: %.10f\n', approx_result);
fprintf('两者误差: %.10e\n', result_original - approx_result);
