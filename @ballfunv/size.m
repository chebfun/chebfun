function varargout = size(v)
% SIZE  Number of expansion coefficients of a BALLFUNV
%   S = SIZE(V) returns the size of the tensor of expansion coefficients 
%   for V, where S = [m,n,p] for an mxnxp tensor of coefficients. 
%   
%   [M, N, P] = SIZE(V) is the same as S = SIZE(V) with M = S(1), N = S(2),
%   and P = S(3). 

% Return the list [m,n,p] corresponding to the size of the ballfun
% functions in f
S = size(v.comp{1});
% Prepare output:
if ( nargout <= 1 ) 
    varargout = { S };
else 
    varargout = { S(1), S(2), S(3) };
end
end
