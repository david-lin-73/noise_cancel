function [est_sig, err_sig, weights] = nlms(noise_ref, mixed_sig, filter_order, step_size)
N = length(mixed_sig);
weights = zeros(filter_order, N);
err_sig = zeros(N, 1);
est_sig = zeros(N, 1);
for i = filter_order:N
    inputs = noise_ref(i:-1:(i-filter_order+1));
    est_sig(i) = weights(:, i)'*inputs;
    err_sig(i) = mixed_sig(i) - weights(:, i)'*inputs;
    %step_size ~ 1e-3
    weights(:, i+1) = weights(:, i) + step_size*e(i)*inputs/(inputs'*inputs);
    %fprintf('The change in w is %d\n', sum((w(:, i+1) - w(:, i)).^2));
end
end