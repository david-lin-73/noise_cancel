function ref = genRectangleRef(length, period, firstPulsePos, pulseWidth)
% Generate rectangle signal with unity amplitude, which is called the reference signal by
% the literature [1]. The output ref is a row vector. The input period is not
% necessarily the interger.
%
% Command:
%    ref = genRectangleRef(length, period, firstPulsePos, pulseWidth)£»
%
% Parameters:
%             ref --- the generated reference signal
%          length --- length of ref
%          period --- the period of pulse, which is not necessarily the interger
%   firstPulsePos --- the first non-zero sample index
%      pulseWidth --- the width of the rectangle width
%
% Author: Zhilin Zhang
% version: 1.0     Date:  Dec.11, 2005


sig = zeros(1,length);
numPeriod = round(length/period)+2;
sig( round([1:numPeriod] * period) ) = 1;
for k = 1 : pulseWidth-1
    sig( round([1:numPeriod] * period)+k ) = 1;
end
sig( 1:round(period)-1 ) = [];

firstPulsePos = round(firstPulsePos);
ref(1: firstPulsePos-1) = zeros(1, firstPulsePos-1);
ref = [ref,sig];
ref(length+1 : end) = [];

