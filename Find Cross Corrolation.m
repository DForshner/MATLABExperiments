function CC = FindCC(Template, Signal)
% FindQRS - Calculates Cross corrolation
% Desc
%
% INPUTS: 
% OUTPUTS: 

% Calcaulate the numerator and denominator of the cross correlation func.
Top = sum( (Template - min(Template)) .* (Signal - min(Signal)) );
Bottom = sum( (Template - min(Template)).^2 ) * sum( (Signal - min(Signal)).^2 );
Bottom = sqrt(Bottom);
CC = Top / Bottom;