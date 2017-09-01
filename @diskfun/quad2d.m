function I = quad2d(f,a,b,c,d)
%QUAD2D  Complete definite integral of DISKFUN. 
%   I = QUAD2D( F ), returns the definite integral of a DISKFUN integrated
%   over its domain of definition. It is the same as sum2(F)
% 
%   I = QUAD2D(F, a, b, c, d), returns the definite integral of a DISKFUN.
%   Integrated over the polar domain [a b] x [c d], where a and b are
%   angles between -pi and pi, and r is a radial parameter between 0 and 1.
% 
% See also INTEGRAL2, SUM2, INTEGRAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Is [a b c d]  subset of the domain of f ? 

%check for empty
if isempty(f)
   I = 0;
   return
end
dom = f.domain;  


% call sum2 if integrating over entire domain
if nargin == 1
    I = sum2(f);
else
if ( ( a < dom(1) ) || ( b > dom(2) ) || ( c < dom(3) ) || ( d > dom(4) ) )
    error('CHEBFUN:DISKFUN:quad2d:domain', ...
        'Can only integrate within the DISKFUN''s domain');
end
f = cart2pol(f, 'cdr'); 
[cols, ~, rows] = cdr(f); %calling cdr after cart2pol makes D = diag(1,1,1,1..1)

%restrict domain
cols = restrict(cols, [c d]);
rows = restrict(rows, [a b]);

% Integrate the rows over their domain.
intRows = sum(rows, [a, b]);

% Create a chebfun of the measure. 
measure = chebfun(@(r) r,[c,d]);

% Multiply the columns by the measure
cols = cols.*(measure*ones(1,size(cols,2)));

% Integrate each column .
intCols = sum(cols, [c d]);

% Put the integrals together to get the final result.
I = sum(intRows.*intCols);


end
