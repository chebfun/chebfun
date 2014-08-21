function varargout = merge(varargin)
%MERGE   Merge information from two CHEBDISCRETIZATION objects.
%   [A, B] = MERGE(A, B) synchronize the properties of CHEBDISCRETIZATIONS A and
%   B so that they will behave compatably when instantiated.
%
%   [A1, A2, ...] = MERGE(A1, A2, ...) is the same as above, but for multiple
%   discretizations.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Initialise:
dom = [];
dimension = 0;
dimAdjust = 0;
projOrder = 0;
varargout = varargin;

% Loop over each discretization:
for k = 1:nargin
    disc = varargin{k};

    dom = domain.merge(dom, disc.domain);
    dimAdjust = mymax(dimAdjust, disc.dimAdjust);
    projOrder = mymax(projOrder, disc.projOrder);

    % Discretizations have to match with the number of subintervals.
    numInt = length(dom) - 1;
    ld1 = length(dimension);
    dimension = [ dimension, repmat(dimension(end), 1, numInt - ld1) ]; %#ok<AGROW>
    ld2 = length(disc.dimension);
    if ( ~isempty(disc.dimension) )
        temp = [ disc.dimension, repmat(disc.dimension(end), 1, numInt - ld2) ];
    end
    dimension = mymax(dimension, temp);
end

for k = 1:nargin
    varargout{k}.domain = dom;
    varargout{k}.dimension = dimension;
    varargout{k}.dimAdjust = dimAdjust;
    varargout{k}.projOrder = projOrder;
end

end

function out = mymax(a, b)
%MYMAX   Take maximum, but avoid empties.
if ( isempty(a) )
    out = b;
elseif ( isempty(b) )
    out = a;
else
    out = max(a, b);
end

end
