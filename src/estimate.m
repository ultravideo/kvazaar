data = dlmread("/dev/stdin", " ");
coeffs = data(1:end, 1:5);
costs = data(1:end, 6);
[beta, sigma, r] = ols(costs, coeffs);
disp(beta)
