function vscl = vscale(f) 
%VSCALE   Vertical scale of a SPHEREFUN.
% 
% VSCL = VSCALE(F) returns the vertial scale of a SPHEREFUN as determined
% by evaluating on a coarse tensor-product grid as defined by SAMPLE
%
% See also SAMPLE.

% NOTE: this is simply a placeholder until vscale@separableApprox has the
% right implementation that does not use chebpolyval2.

% TODO: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% If f is an empty SPHEREFUN, VSCL = 0: 
if ( isempty( f ) ) 
    vscl = 0; 
    return
end

% Get the degree of the SPHEREFUN:
[m, n] = length(f); 

% If F is of low degree, then oversample: 
m = min(max(m, 9),2000); 
n = min(max(n, 9),2000); % cannot afford to go over 2000x2000. 

% Calculate values on a tensor grid: 
vals = sample(f, m, n); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end