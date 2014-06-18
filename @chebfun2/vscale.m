function vscl = vscale(f) 
%VSCALE   Vertical scale of a CHEBFUN2.
% 
% VSCL = VSCALE(F) returns the vertial scale of a CHEBFUN2 as determined
% by evaluating on a coarse Chebyshev tensor-product grid. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% If f is an empty Chebfun2, VSCL = 0: 
if ( isempty(f) ) 
    vscl = 0; 
    return
end

% Get the degree of the CHEBFUN2:
[m, n] = length(f); 

% If F is of low degree, then oversample: 
m = max(m, 9); 
n = max(n, 9); 

% Calculate values on a tensor grid: 
vals = chebpolyval2(f, m, n); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end