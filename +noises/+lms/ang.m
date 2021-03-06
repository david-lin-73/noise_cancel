function [est_sig, err_sig, weights] = ang(ref_sig, mixed_sig, filt_ord, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description:
%
% This is an implementation of the adaptive step size LMS derived from
% Wee-Peng Ang's work. This is the simplified Beneveniste's common step
% size (Linear) 
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=912925
%
% Arguments (required):
%
% ref_sig       [Nx1]       Reference signal of length N  
% mixed_sig     [Nx1]       Mixed signal (desired) of length N
% filt_ord                  Number of desired filter taps
%               
% Arguements (optional):
% rho                       Estimate of gradient power (should be small)
% step_size                 Initial learning rate (should be small)
% plot                      plots the est, err, mixed, and reference sigs
%
% Outputs
% est_sig       [Nx1]       Estimated signal of reference through filter
% err_sig       [Nx1]       Error from mixed_signal - est_signal    
% weights       [LxN]       L tap weights (filter order) over entire signal
%                                of length N
%
% Notes:
% Small filter size (<10) with extremely small rho (<10^-7) works well
% Using noise as input does not work as well as actual estimate of what the
% mixed signal should end up producing
% rho:  chosen such that the excess Mean Square Error (MSE) and temporal
% rms excess (MSE) are the same larger rho = faster convergence and lower
% excess MSE but produces large fluctuation of temporal MSE around avg val
%
% Future work:
% Implementing a limitation to mu that is hard limited to a certain value
% to ensure 
% Can also do multiple step size for each tap as well as multiplicative
% Work History
% 2/5/2020          Initial comments and base code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% extract optional UI 
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'step_size'; step_size = varargin{arg+1};
            case 'rho'; rho = varargin{arg+1};
            case 'display'; display = 1;
        end
    end
end


% set default parameters otherwise
if ~exist('rho', 'var'); rho = 4e-9; end
if ~exist('step_size', 'var'); step_size = 1/((ref_sig'*ref_sig)); end
if ~exist('display', 'var'); display = 0; end


% empirical constant 0 < alpha < 1
alpha = 0.95;

% initialize various parameters
N = length(mixed_sig);
est_sig = zeros(N, 1);
weights = zeros(filt_ord, N);
err_sig = zeros(N, 1);
steps = zeros(N, 1);
phi = zeros(filt_ord, N);
input_track = zeros(filt_ord, N);

% set first step_size for algorithm
steps(filt_ord-1) = step_size;

% start algorithm
for i = filt_ord:N-1
    % take a flipped subset of the reference 
    %inputs = flipud(ref_sig((i-filt_ord+1):i));
    inputs = ref_sig(i:-1:(i-filt_ord+1));
    input_track(:, i) = inputs;
    
    % calculate the estimated and error signal respectively
    est_sig(i) = weights(:, i)'*inputs;
    err_sig(i) = mixed_sig(i) - weights(:, i)'*inputs;
    
    %%% Ang's alg: change the step size
    phi(:, i) = alpha*phi(:, i-1) + err_sig(i-1)*input_track(:, i-1);
    steps(i) = steps(i-1) + rho*err_sig(i)*inputs'*phi(:, i);  
    
    % update the weights
    weights(:, i+1) = weights(:, i) + steps(i)*err_sig(i)*inputs;

end
if display == 1
    figure
    t = (1:N);
    subplot(4,1,1)
    plot(est_sig)
    title('Estimated signal')
    subplot(4,1,2)
    plot(t, err_sig)
    title('Error Signal')
    subplot(4,1,3)
    plot(t, ref_sig)
    title('Reference signal')
    subplot(4, 1, 4)
    plot(t, mixed_sig)
    title('Mixed Signal')
end
end


