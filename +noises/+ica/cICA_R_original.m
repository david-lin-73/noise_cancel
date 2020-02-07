function [Y, w] = cICA_R(X, R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3984144/pdf/pone.0094211.pdf

% CHECK MEANS
% assuming X is a whitened matrix that is MxN where M is the number of RV
% and N is the number of samples
[M, N] = size(X);

% set parameters 
gam = 1;
thresh = 0.9;
lambda = 1;
mu = 1;
learnrate = 1;
EndCond = 1e-9;
MaxIter = 400;

a1 = 1; % but a1 exists in the sapce [1, 2]

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
%CorrMat = X(:,[1:end]) * X(:,[1:end])' / N;
flag = 1;
iter = 0;
wold = w;
while flag == 1
    Y = w*X;
    % gaussian variable with same parameters as Y zero_mean unit variance
    V = normrnd(mean(Y), std(Y), size(Y));
    % can add a + weight constant in front to make certain contrast calcs
    %   more significant (w_c) should be a vector for matrix W 
    w_c = 1;
    %sign(rho) is necessary IF changing threshold 
    rho = (w_c * mean(Neg(Y) - Neg(V)));
    wJ = rho * X*Neg_1(Y)'/N;

    % not expected value
    %wG = 1/2 * mean(X*mean(Close_1(Y, R)), 2);s
    wG = 1/2 * mu * X*(Close_1(Y, R))'/N; 
   
    wH = lambda * (X*Y')/N;
    
    %1st derivative of Lagrangian
    wL = wJ - wG - wH;
    
    % might be a summation check tIME sPAce corr paper
    % 2nd derivative of Lagrangian 
    w2L = rho*mean(Neg_2(Y)) - mu/2 * mean(Close_2(Y, R)) - lambda;
   
    % update components
    
    w = (w' - inv(CorrMat)*(learnrate * wL / w2L))';
    w = w/norm(w);
    % update 1/26 added a changing threshold improves convergence
    t = thresh*log(exp(1)+iter);
    %mu = max(0, mu + gam * (mean(Close(Y, R)) - thresh));
    mu = max(0, mu + gam * (mean(Close(Y, R)) - t));
    lambda = lambda + gam * (mean(Y.^2) - 1);
    %%% TO DO: LOOK AT HOW TO DECORRELATE %%%
    %w = (w*w').^(-1/2)*w;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dot product = 1 (means that they point in the same direction
    deltw = abs(abs(w*wold') - 1);
    fprintf('at %d iterw changed by: %g\n',iter, deltw);
    if (abs(abs(w*wold') - 1)) < EndCond && (mean(Close(Y, R)) - t) <= 0
    %if (abs(abs(w*wold') - 1)) < EndCond && (mean(Close(Y, R)) - thresh) <= 0
        fprintf('at iter %n, w converged to: %g\n',iter, w);
        flag = 0;
    end
    
    if iter > MaxIter
        flag = 0;
    end
    iter = iter + 1;
    wold = w;
end
end

%z = w*X;
%plot(z(201:700)); axis([-inf,inf,-4,4]);

    
    
   
    
    
    
    
    
    


