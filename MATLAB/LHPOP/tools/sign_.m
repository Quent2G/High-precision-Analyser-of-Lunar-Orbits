%--------------------------------------------------------------------------
% sign: returns absolute value of a with sign of b
%
% Last modified:   2018/01/27   M. Mahooti
%--------------------------------------------------------------------------
function [result] = sign_(a, b)

if (b>=0)
    result = abs(a);
else
    result = - abs(a);
end