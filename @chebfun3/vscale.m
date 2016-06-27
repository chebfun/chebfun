function vscl = vscale(f) 
%VSCALE   Vertical scale of a CHEBFUN3. 
% VSCL = VSCALE(F) returns the vertial scale of a CHEBFUN3 object F as 
%   determined by evaluating F on a coarse tensor-product grid. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% [TODO]: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% If f is an empty CHEBFUN3, VSCL = 0: 
if ( isempty(f) ) 
    vscl = 0; 
    return
end

techCol = get(f.cols.funs{1}, 'tech');

% Get the degree of the CHEBFUN3:
[m, n, p] = length(f);

% If F is of low degree, then oversample: 
m = min(max(m, 9), 41); 
n = min(max(n, 9), 41); 
p = min(max(p, 9), 41); % cannot afford to go over 41x41x41. 

% Calculate values on a tensor grid: 
dom = f.domain;
x = mypoints(m, dom(1:2), techCol);
y = mypoints(n, dom(3:4), techCol);
z = mypoints(p, dom(5:6), techCol);
[xx, yy, zz] = ndgrid(x, y, z);
vals = feval(f, xx, yy, zz); 

% Take the absolute maximum: 
vscl = max(abs(vals(:)));

end

%%
function x = mypoints(n, dom, tech)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
if ( isa(tech(), 'chebtech2') )    
    x = chebpts(n, dom, 2);
elseif ( isa(tech(), 'chebtech1') )    
    x = chebpts(n, dom, 1);
elseif ( isa(tech(), 'trigtech') )
    x = trigpts(n, dom);
else
    error('CHEBFUN:CHEBFUN3:vscale:mypoints:techType', ...
        'Unrecognized technology');
end

end