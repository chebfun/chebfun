function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] horizontally concatenates the TRIGTECH objects A and B to form an
%   array-valued TRIGTECH. [A,B] does the same. Any number of TRIGTECH objects
%   can be concatenated within one pair of brackets. Vertical concatenation is
%   not supported.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Remove empties:
empties = cellfun(@isempty, varargin);
if ( all(empties) )
    out = varargin{1};
    return
else
    varargin(empties) = [];
end

% Check that all objects are indeed TRIGTECHs.  If not then we should not
% concatenate.
if ~all( cell2mat(cellfun(@(f) isa(f,'trigtech'), varargin, 'UniformOutput', false)) )
    error('CHEBFUN:TRIGTECH:horzcat:typeMismatch','Incompatible operation between objects. Make sure functions are of the same type.');
end

% Prolong each TRIGTECH to the same length:
n = max(cellfun(@length, varargin));
F = cellfun(@(f) prolong(f, n), varargin, 'UniformOutput', false);

% Extract the data and collate the an array-valued TRIGTECH:
out = varargin{1};

% Coeffs and Values:
out.coeffs = cell2mat(cellfun(@(f) f.coeffs, F, 'UniformOutput', false));
out.values = cell2mat(cellfun(@(f) f.values, F, 'UniformOutput', false));

% Vscales:
vscales = cellfun(@(f) f.vscale, F, 'UniformOutput', false);
out.vscale = cell2mat(vscales);

% Epslevel:
epslevels = cellfun(@(f) f.epslevel, F, 'UniformOutput', false);
out.epslevel = cell2mat(epslevels);

% Hscale:
out.hscale = max(cellfun(@(f) f.hscale, F));

% IsReal:
areReal = cellfun(@(f) f.isReal, F, 'UniformOutput', false);
out.isReal = cell2mat(areReal);

end
