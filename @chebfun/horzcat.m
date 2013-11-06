function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN objects.
%   [A B] is the horizontal concatenation of CHEBFUN objects A and B. [A,B] is
%   the same thing. Any number of CHEBFUN objects can be concatenated within
%   one pair of brackets. Vertical concatenation is not supported.
%
% See also VERTCAT, CAT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Vertical concatenation.

% Promote doubles to CHEBFUN objects:
chebfunLocs = cellfun('isclass', varargin, 'chebfun');
dom = varargin{find(chebfunLocs, 1, 'first')}.domain;
doubleLocs = find(~chebfunLocs);
for k = doubleLocs
    varargin{k} = chebfun(varargin{k}, dom);
end

% Ensure that the domains match:
if ( any(diff(cellfun(@(f) length(f.domain), varargin))) )
    error('CHEBFUN:horzcat:noSupport', ...
        'HORZCAT of a CHEBFUNs on different domains is not yet supported.');
end
tol = max(cellfun(@(f) hscale(f).*epslevel(f), varargin));

if ( any(cellfun(@(f) any(f.domain - dom) > tol, varargin)) )
    error('CHEBFUN:horzcat:noSupport', ...
        'HORZCAT of a CHEBFUNs on different domains is not yet supported.');
end

% Concatenate the FUNs:
out = varargin{1};
for k = 1:numel(dom)-1
   funs = cellfun(@(f) f.funs{k}, varargin, 'UniformOutput', false);
   out.funs{k} = horzcat(funs{:});
end

end
