function r = roots( F, varargin )
%ROOTS   Find the common zeros of a DISKFUNV object.
%   r = ROOTS(F) finds the common zeros of the two bivariate functions F(1) and
%   F(2) in their domain of definition under the assumption that the solution
%   set is zero-dimensional. R is a matrix with two columns storing the x- and
%   y-values of the solutions. 
%
% See also  CHEBFUN2V/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


%call roots in diskfun
if isempty(F)
    r = [];
    return
end

r = roots(F.components{1}, F.components{2}, varargin{:}); 