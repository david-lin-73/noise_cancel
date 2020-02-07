function Weights = fpfica(data, num_comp) 
% Description: 
%                   Fixed point ICA 
%
% Inputs:   
%             data      - NxM dataset where N is the number of observables
%                           (i.e random variables) and M is the number of 
%                           samples per variable. DATA IS ASSUMED TO BE 
%                           ZERO MEAN AND WHITENED see preprocessing.m
%             num_comp  - Number of components that you wish to
%                            extract from the data
%
% Outputs:          
%             Weights   - MxC weight matrix where M is the number of
%                           samples per observable and C is the number of 
%                           components to extract. 
%                           Eventually Weights'*data = components
%                   
%
%
%


% contrast functions based on negentrophy
G = @(w) tanh(w);
G_ = @(w) 1 - G(w).^2;
%%% more contrast functions if desired G is the function G_ is its
%%% derivative
% G = @(Sk) Sk .* exp(-0.5 * Sk.^2);
% G_ = @(Sk) (1 - Sk.^2) .* exp(-0.5 * Sk.^2);
%G = @(w) -1 * exp(-0.5*(w.^2));
%G_ = @(w) w.*exp(-0.5*(w.^2))

% 
W = rand(size(data, 1), num_comp);
for i = 1:num_comp
    w = W(:, i);
    w = w/norm(w);
    iter = 0; 
    thresh = inf;
    while iter < 100 && thresh > 1e-10

        wnew = mean(data*G(w'*data)', 2) - mean(G_(w'*data)) *w;


        wdecor = zeros(size(wnew));
        for j = 1:i-1
            wdecor = wdecor + wnew'*W(:, j)*W(:, j);
        end
        wnew = wnew - wdecor;

        wnew = wnew/norm(wnew);
        iter = iter + 1;
        thresh = abs(abs(w'*wnew) - 1);
        w = wnew;
    end
    W(:, i) = wnew;
end
Weights = W
end