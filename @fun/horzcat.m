function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] horizontally concatenates the FUN objects A and B to form an
%   array-valued FUN. [A,B] does the same. Any number of FUN objects can be
%   concatenated within one pair of brackets. Vertical concatenation is not
%   supported.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Remove empties:
empties = cellfun(@isempty, varargin);
if ( all(empties) )
    out = varargin{1};
    return
else
    varargin(empties) = [];
end

% Make an output FUN:
out = varargin{1};

% Extract the ONEFUNs:
onefuns = cellfun(@(f) f.onefun, varargin, 'UniformOutput', false);

% Concatenate the ONEFUNs:
out.onefun = horzcat(onefuns{:});

end
