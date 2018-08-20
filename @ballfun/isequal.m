function b = isequal(f, g)
% ISEQUAL Test the equality between two BALLFUN functions
%   ISEQUAL(f, g) is the boolean f == g

% Test if f = g
if (nnz(size(f)-size(g))==0)
    b = iszero(f-g);
else
    error('BALLFUN:isequal:unknown', ...
    ['Undefined function ''isequal'' for different size of ballfun functions : ' ...
     '%s and %s.'], mat2str(size(f)), mat2str(size(g)));
end
end
