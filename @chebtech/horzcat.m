function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] horizontally concatenates the CHEBTECH objects A and B to form an
%   array-valued CHEBTECH. [A,B] does the same. Any number of CHEBTECH objects
%   can be concatenated within one pair of brackets. Vertical concatenation is
%   not supported.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Remove empties:
empties = cellfun(@isempty, varargin);
if ( all(empties) )
    out = varargin{1};
    return
else
    varargin(empties) = [];
end

% Check that all objects are CHEBTECHs. If not we cannot concatenate.
if ( ~all( cellfun(@(f) isa(f, 'chebtech'), varargin) ) )
    error('CHEBFUN:CHEBTECH:horzcat:typeMismatch', ...
    'Incompatible concatenation. Ensure discretizations are of the same type.');
end

% Prolong each Chebtech to the same length:
n = max(cellfun(@length, varargin));
F = cellfun(@(f) prolong(f, n), varargin, 'UniformOutput', false);

% Extract the data and collate the an array-valued CHEBTECH:
out = varargin{1};

% Coeffs:
out.coeffs = cell2mat(cellfun(@(f) f.coeffs, F, 'UniformOutput', false));

end
