function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN objects.
%   [A B] horizontally concatenates the CHEBFUN objects A and B to form an
%   array-valued CHEBFUN. [A,B] does the same. Any number of CHEBFUN objects
%   can be concatenated within one pair of brackets. Vertical concatenation is
%   not supported.
%
% See also VERTCAT, CAT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Vertical concatenation. (Chebmatrix / quasimatrix).

% [TODO]: Currently if breakpoints don't match then we restrict. In future this
% should instead return a quasimatrix, and the documentation above should be
% updated to reflect this.

% Remove empties:
empties = cellfun(@isempty, varargin);
if ( all(empties) )
    out = varargin{1};
    return
else
    varargin(empties) = [];
end

% Promote doubles to CHEBFUN objects:
chebfunLocs = cellfun('isclass', varargin, 'chebfun');
domain1 = varargin{find(chebfunLocs, 1, 'first')}.domain;
doubleLocs = find(~chebfunLocs);
for k = doubleLocs
    varargin{k} = chebfun(varargin{k}, domain1);
end

% Grab the domains of each of the inputs:
allDomainsCell = cellfun(@(f) f.domain, varargin, 'UniformOutput', false);

% Ensure that the domains match:
domainEnds = allDomainsCell{1}([1 end]);
if ( any(cellfun(@(d) any(d([1 end]) ~= domainEnds), allDomainsCell)) )
    error('CHEBFUN:horzcat:domains', 'Inconsistent domains.');
end

% Check to see if interior breakpoints differ:
differentBreakpoints = false;
if ( any(diff(cellfun(@(d) length(d), allDomainsCell))) )
    differentBreakpoints = true;
else
    tol = max(cellfun(@(f) hscale(f).*epslevel(f), varargin));
    if ( any(cellfun(@(d) any(d - allDomainsCell{1}) > tol, allDomainsCell)) )
        differentBreakpoints = true;
    end
end
% Restrict / overlap if they do:
if ( differentBreakpoints )
    unionOfDomains = unique([allDomainsCell{:}]);
    varargin = cellfun(@(f) restrict(f, unionOfDomains), varargin, ...
        'UniformOutput', false);
end

% Concatenate the FUNs:
out = varargin{1};
numInts = numel(out.domain) - 1;
for k = 1:numInts
   funs = cellfun(@(f) f.funs{k}, varargin, 'UniformOutput', false);
   out.funs{k} = horzcat(funs{:});
end

end
