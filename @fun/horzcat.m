function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] is the horizontal concatenation of FUN objects A and B on the same
%   domain. [A,B] is the same thing. Any number of FUN objects can be
%   concatenated within one pair of brackets. Vertical concatenation is not
%   supported.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Make an output FUN:
out = varargin{1};

tol = get(out, 'hscale')*get(out, 'epslevel');
if ( any( cellfun(@(f) any(f.domain - out.domain) > tol, varargin)) )
    error('FUN:horzcat:noSupport', 'Domain mismatch.');
end

% Extract the ONEFUNs:
onefuns = cellfun(@(f) f.onefun, varargin, 'UniformOutput', false);

% Concatenate the ONEFUNs:
out.onefun = horzcat(onefuns{:});

end