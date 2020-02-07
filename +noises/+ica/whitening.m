function [Wtrans, Wmatrix] = whitening(CentMatrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Description:
%
% This function whitens the inputed matrix which takes the input matrix
% whose rows are different random variables (RV) and transforms it into a
% new set of RV's such that the covariance is now the identity matrix
% meaning: Uncorrelated and unit variance. Used in preprocessing before ICA
% 
% Input (required)
% 
% CentMatrix    [NxM]    N random variables of M samples which have
%                           already been zero-meaned
%
% Ouputs:
%
% Wtrans        [NxM]    N uncorreleated random variables of M samples
%                           which now have unit variance
% Wmatrix       [NxN]    whitening matrix that converts input into whitened
%                           matrix
% Work History
% 2/5/2020          Initial comments and base code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find covariance of input
CovMatrix = cov(CentMatrix');

% Calculate whitening matrix
[U, S, ~] = svd(CovMatrix, 'econ');
D = diag(diag(1./sqrt(S)));
Wmatrix = U*D*U';

% Apply the whitening matrix on input
Wtrans = Wmatrix*CentMatrix;
end