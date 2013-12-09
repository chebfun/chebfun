function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] horizontally concatenates the CHEBTECH objects A and B to form an
%   array-valued CHEBTECH. [A,B] does the same. Any number of CHEBTECH objects
%   can be concatenated within one pair of brackets. Vertical concatenation is
%   not supported.

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
  
% Prolong each Chebtech to the same length:
n = max(cellfun(@length, varargin));
F = cellfun(@(f) prolong(f, n), varargin, 'UniformOutput', false);

% Extract the data and collate the an array-valued CHEBTECH:
out = varargin{1};

% Coeffs and Values:
out.coeffs = cell2mat(cellfun(@(f) f.coeffs, F, 'UniformOutput', false));
out.values = cell2mat(cellfun(@(f) f.values, F, 'UniformOutput', false));

% Vscales:
vscales = cellfun(@(f) f.vscale, F, 'UniformOutput', false);
out.vscale = max(cell2mat(vscales));
vscales = cellfun(@max, vscales);
%vscales = max(max(vscales)); 

% Epslevel:
out.epslevel = max(cellfun(@(f) f.epslevel, F).*vscales)./max(vscales);

% Hscale:
out.hscale = max(cellfun(@(f) f.hscale, F));

end
