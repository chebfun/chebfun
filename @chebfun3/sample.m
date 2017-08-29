function varargout = sample(f, varargin)
%SAMPLE   Values of a CHEBFUN3 object on a tensor product grid.
%   X = SAMPLE(F) returns the tensor of values of F on a tensor product 
%   grid.
%
%   [CORE, C, R, T] = SAMPLE(F) returns the low rank representation of the
%   values of F on a tensor product grid so that X = CORE x_1 C x_2 R x_3 T.
%
%   [CORE, C, R, T] = SAMPLE(F, M, N, P) returns the values of F on an
%   M x N x P tensor product grid.
%
% See also CHEBFUN3/FEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(f) )
    varargout = {[]}; 
    return
end

if ( nargin == 4 )
     m = varargin{1};
     n = varargin{2};
     p = varargin{3};
else
    [m,n,p] = length(f);
    m = max(m, 51);
    n = max(n, 51);
    p = max(p, 51);
end

% Use Slice-Tucker decomposition so we can keep it in low rank form:
[fCore, fCols, fRows, fTubes] = tucker(f);
Cvals = sample(fCols, m);
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