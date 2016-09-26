function vscl = vertscale(f)
%VSCALE   Vertical scale of a CHEBFUN3T object.
%   VSCL = VSCALE(F) returns the vertial scale of a CHEBFUN3T as determined
%   by evaluating on a coarse tensor-product grid. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If f is an empty CHEBFUN3T, VSCL = 0: 
if ( isempty(f) ) 
    vscl = 0; 
    return
end

% Get the degree of the CHEBFUN3T:
[m, n, p] = length(f);

% If F is of low degree, then oversample: 
m = min(max(m, 9), 41); 
n = min(max(n, 9), 41); 
p = min(max(p, 9), 41); % cannot afford to go over 41 x 41 x 41. 

% Calculate values on a tensor grid: 
dom = f.domain;

x = chebpts(m, dom(1:2));
y = chebpts(n, dom(3:4));
z = chebpts(p, dom(5:6));
[xx, yy, zz] = ndgrid(x, y, z);
vals = feval(f, xx, yy, zz); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end