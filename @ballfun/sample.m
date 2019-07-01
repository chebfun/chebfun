function X = sample( f, varargin )
%SAMPLE      Samples f on a tensor product grid.
%   X = SAMPLE(F) returns the matrix of values of F on an m-by-n-by-p tensor
%   product grid, where m is length of columns of F, n is length of the 
%   rows, and p is the length of tubes.
%
%   X = SAMPLE(F, M, N, P) returns the matrix of values of F on an n-by-m-by-p 
%   tensor product grid.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(f) )
    varargout = {[]}; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n, p] = size(f);
elseif ( nargin == 2 || nargin == 3 ) 
    error('CHEBFUN:BALLFUN:sample:inputs', 'Dimension not specified.'); 
else
    m = varargin{1}; 
    n = varargin{2};
    p = varargin{3};
    if ( (m <= 0) || (n <= 0) || (p <= 0) )
        error('CHEBFUN:BALLFUN:sample:inputs', ['Number of sample ' ...
             'points must be positive.']);
    end
end

% Evaluation points
r = chebpts(2*m-1);
r = r(m:end);
l  = pi*trigpts(n);
t = pi*linspace(0,1,p);

% Evaluate f
X = fevalm(f, r, l, t);
end