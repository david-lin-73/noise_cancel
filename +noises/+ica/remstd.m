function sig = remstd(OldSig)
% remstd --- make the input signals' variance unit
% The input OldSig is in the form N*P, where N is the number of signals, and P is
% the length of each signal, then the returned sig is also in the form N*P,
% but each row (each signal) has unit variance.
%
% Command:
%    sit = remstd( OldSig );
%
% See also:
%     centering    whiten    standarize
%
% Author: Zhilin Zhang
%
% Last version: 1.5     Date: Apr.17,2005
%      Version: 1.0     Date: Oct.31,2003


n=size(OldSig,1);
for t=1:n
    sig(t,:)=OldSig(t,:)/std(OldSig(t,:));
end
