function [est_sig, err_sig, weights] = lms(noise_ref, mixed_sig, filter_order, step_size)
N = length(mixed_sig);
weights = zeros(filter_order, N);
err_sig = zeros(N, 1);
est_sig = zeros(N, 1);
%step_size = 1/(2*(noise_ref'*noise_ref))
for i = filter_order:N
    inputs = noise_ref(i:-1:(i-filter_order+1));
    est_sig(i) = weights(:, i)'*inputs;
    err_sig(i) = mixed_sig(i) - weights(:, i)'*inputs;

    weights(:, i+1) = weights(:, i) + step_size*err_sig(i)*inputs;
    fprintf('The change in w is %d\n', sum((weights(:, i+1) - weights(:, i)).^2));
end
end