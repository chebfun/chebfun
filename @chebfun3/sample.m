function varargout = sample(f, varargin)
%SAMPLE     Values of f on a tensor product grid.
%
%   X = SAMPLE(F) returns the tensor of values of F on a tensor product 
%   grid.
%
%   [CORE, C, R, T] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = CORE x_1 C x_2 R x_3 T.
%
%   [CORE, C, R, T] = SAMPLE(F,M,N,P) returns the values of F on a
%   M-by-N-by-P tensor product grid.

% Empty check. 
if ( isempty(f) )
    varargout = {[]}; 
    return
end

[m,n,p] = length(f);
m = max(m, 50);
n = max(n, 50);
p = max(p, 50);

% For now just call chebpolyval3.  In the future, this function should call
% an appropriate sample function on the underlying tech type. This would
% require adding "sample" functions at the tech level.
%[fCore, fCols, fRows, fTubes] = chebpolyval3(f, varargin{:});


% Use ST decomposition so we can keep it in low rank form: 
[fCore, fCols, fRows, fTubes] = st(f);
Cvals = sample(fCols, m);
%Cvals = sample(chebfun(C,'splitting', 'on'), m);
Rvals = sample(fRows, n);
Tvals = sample(fTubes, p);

% Evaluate: 
if ( nargout <= 1 )
    varargout = {chebfun3.txm(chebfun3.txm(chebfun3.txm(fCore, Cvals, ...
        1), Rvals, 2), Tvals, 3)};
else
    varargout = {fCore, Cvals, Rvals, Tvals};
end

end