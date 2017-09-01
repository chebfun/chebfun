function vscl = vscale(f) 
%VSCALE   Vertical scale of a SEPARABLEAPPROX.
% 
% VSCL = VSCALE(F) returns the vertial scale of a SEPARABLEAPPROX as determined
% by evaluating on a coarse tensor-product grid. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% If f is an empty SEPARABLEAPPROX, VSCL = 0: 
if ( isempty( f ) ) 
    vscl = 0; 
    return
end

% Get the degree of the SEPARABLEAPPROX:
[m, n] = length(f); 

% If F is of low degree, then oversample: 
m = min(max(m, 9),2000); 
n = min(max(n, 9),2000); % cannot afford to go over 2000x2000. 

% Calculate values on a tensor grid: 
vals = sample(f, m, n); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end
