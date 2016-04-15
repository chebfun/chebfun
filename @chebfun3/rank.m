function varargout = rank(F)
%RANK   Multilinear rank of a CHEBFUN3 object F.
%
%   If three outputs are asked for, then it is the size of core tensor, 
%   i.e., r = [size(F.cols,2), size(F.rows,2), size(F.tubes,2)];
%   If just one output is asked for, then r is the max size of the core 
%   tensor. 

r = size(F.core);
if numel(r)<3
    r = [r, 1];
end

% Output:
if ( nargout <= 1 )
    varargout = {max(r(:))};
else
    varargout = {r(1), r(2), r(3)};
end

end