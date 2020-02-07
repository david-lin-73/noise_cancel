function [Y, w] = cICA_R(X, R, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description
% cICA based off of Negentropy constrained optimization to extract desired
% Independent Component (IC) from mixture of observations. It uses a
% steepest descent algorithm by converting the optimization equation into a
% Lagrangian equation and solving until convergence. References used are
% listed below: 
% (Wang)
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3984144/pdf/pone.0094211.pdf
% (Lu)
% http://read.pudn.com/downloads133/doc/project/567667/ICA_with_Reference.pdf
% (Zhang)
% http://dsp.ucsd.edu/~zhilin/papers/McICA.pdf
%
% Inputs (required)
% X         [NxM]       N mixed signals of M samples each
% R         [1xM]       Reference signal of M samples (1 signal for now)
% 
% Inputs (optional)
% sign                  Determines whether to use sign(rho) in loop
% gamma                 

% Notes: 
% The closeness measure might require the output and reference be
% normalized (Lu 2249)
% Future work:
% Adding capability for multiple references
% Adding multiple weights (decorrelation at each iteration) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assuming X is a whitened matrix that is MxN where M is the number of RV
% and N is the number of samples
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'sign'; Sign = true;
        end
    end
end

if ~exist('sign', 'var');       Sign = false;       end
if ~exist('gamma', 'var');      gam = 1;            end 
if ~exist('threshold', 'var');  thresh = 0.9;       end
if ~exist('lambda', 'var');     lambda = 1;         end 
if ~exist('learnrate', 'var');  learnrate = 0.98;   end
if ~exist('converge', 'var');   EndCond = 1e-9;     end 
if ~exist('iteration', 'var');  MaxIter = 400;      end 

[M, N] = size(X);



% set parameters 
gam = 1;
thresh = 0.9;
lambda = 1;
mu = 1;
learnrate = 0.98;
EndCond = 1e-9;
MaxIter = 400;

a1 = 1; % but a1 exists in the space [1, 2]

% The negentropy function 
Neg = @(m) log(cosh(a1*m));
Neg_1 = @(m) tanh(m);
Neg_2 = @(m) 1 - tanh(m).^2;

% the correlation function (Zhang)
Close = @(o, r) (o-r).^2;
Close_1 = @(o, r) 2*(o-r);
Close_2 = @(o, r) 2;

% column row vec, eventually will implement complete weight matrix
% W will be a MxC matrix where M is the numbmer of RV and 
% C is the number of components to extract st Y = WX; 
w = ones(1, size(X, 1));
w = w/norm(w);
% dummy for now eventually this will be our reference signal
%R = ones(size(w*X));
CorrMat = cov(X');


flag = 1;
iter = 0;
wold = w;
while flag == 1
    Y = w*X;
    % gaussian variable with same parameters as Y zero_mean unit variance
    V = normrnd(mean(Y), std(Y), size(Y));
    % can add a + weight constant in front to make certain contrast calcs
    %   more significant (w_c) should be a vector for matrix W (
    w_c = 1;
    %sign(rho) is necessary IF changing threshold 
    rho = (w_c * mean(Neg(Y) - Neg(V)));
    if Sign == true
        rho = sign(rho);
    end
    
    % partial derivative of negentropy
    wJ = rho * X*Neg_1(Y)'/N;

    % partial derivative of closeness constraint
    wG = 1/2 * mu * X*(Close_1(Y, R))'/N; 
   
    % partial derivative of 
    wH = lambda * (X*Y')/N;
    
    %1st derivative of Lagrangian
    wL = wJ - wG - wH;
    
    % might be a summation check tIME sPAce corr paper
    % 2nd derivative of Lagrangian 
    w2L = rho*mean(Neg_2(Y)) - mu/2 * mean(Close_2(Y, R)) - lambda;
   
    % update components
    % Update 2/5/2020 variable learn rate can be added: learnrate^i (Wang)
    w = (w' - inv(CorrMat)*(learnrate * wL / w2L))';
    w = w/norm(w);
    
    % start with small thresh then increase to globally converge (Lu)
    % update 1/26 added a changing threshold improves convergence
    t = thresh*log(exp(1)+iter);

    mu = max(0, mu + gam * (mean(Close(Y, R)) - t));
    lambda = lambda + gam * (mean(Y.^2) - 1);

    % dot product = 1 (means that they point in the same direction
    deltw = abs(abs(w*wold') - 1);
    fprintf('at %d iterw changed by: %g\n',iter, deltw);
    
    % check if new weights have converged and output correlates with ref
    if (deltw) < EndCond && (mean(Close(Y, R)) - t) <= 0
        fprintf('at iter %n, w converged to: %g\n',iter, w);
        flag = 0;
    end
    
    % check if maximum convergence 
    if iter > MaxIter
        flag = 0;
        fprintf('The algorithm has stopped reaching %d iterations', MaxIter')
    end
    
    %update for next iteration
    iter = iter + 1;
    wold = w;
end
end



    
    
   
    
    
    
    
    
    


