function vscl = vscale(f) 
%VSCALE   Vertical scale of a CHEBFUN3.
% 
% VSCL = VSCALE(F) returns the vertial scale of a CHEBFUN3 as determined
% by evaluating on a coarse tensor-product grid. 

% TODO: Should this also be taking the maximum along the edges when we are
% evaluating at 1st kind grids. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If f is an empty CHEBFUN3, VSCL = 0: 
if ( isempty(f) ) 
    vscl = 0; 
    return
end

% Check underlying tech is a TRIGTECH: 
techCol = get(f.cols.funs{1}, 'tech');
techRow = get(f.rows.funs{1}, 'tech');
techTube = get(f.tubes.funs{1}, 'tech');
if ( isa(techCol(), 'trigtech') && isa(techRow(), 'trigtech') ...
        && isa(techTube(), 'trigtech') )
    tech = trigtech;
else
    tech = chebtech2;
end


% % Get the degree of the CHEBFUN3:
[m, n, p] = length(f);

% If F is of low degree, then oversample: 
m = min(max(m, 9), 41); 
n = min(max(n, 9), 41); 
p = min(max(p, 9), 41); % cannot afford to go over 40x40x40. 

% Calculate values on a tensor grid: 
%vals = sample(f, m, n, p); 
dom = f.domain;
x = mypoints(m, dom(1:2), tech);
y = mypoints(n, dom(3:4), tech);
z = mypoints(p, dom(5:6), tech);
[xx, yy, zz] = ndgrid(x, y, z);
% [xx, yy, zz] = ndgrid(chebpts(m, dom(1:2)), chebpts(n, dom(3:4)), ...
%     chebpts(p, dom(5:6)));
vals = feval(f, xx, yy, zz); 

% Take the absolute maximum: 
vscl = max(abs(vals(:))); 

end

%%
function x = mypoints(n, dom, tech)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
%tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts( n, dom, 2 ); 
elseif ( isa(tech, 'chebtech1') )
    x = chebpts( n, dom, 1 ); 
elseif ( isa(tech, 'trigtech') )
    x = trigpts( n, dom ); 
else
    error('CHEBFUN:CHEBFUN3:vscale:mypoints:techType', ...
        'Unrecognized technology');
end

end