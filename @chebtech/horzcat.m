function out = horzcat(varargin)
%HORZCAT   Horizontal concatenation.
%   [A B] horizontally concatenates the CHEBTECH objects A and B to form an
%   array-valued CHEBTECH. [A,B] does the same. Any number of CHEBTECH objects
%   can be concatenated within one pair of brackets. Vertical concatenation is
%   not supported.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
  
% Prolong each Chebtech to the same length:
n = max(cellfun(@length, varargin));
F = cellfun(@(f) prolong(f, n), varargin, 'UniformOutput', false);

% Extract the data and collate the an array-valued CHEBTECH:
out = varargin{1};
out.coeffs = cell2mat(cellfun(@(f) f.coeffs, F, 'UniformOutput', false));
out.values = cell2mat(cellfun(@(f) f.values, F, 'UniformOutput', false));
out.vscale = cellfun(@(f) f.vscale, F);
out.hscale = max(cellfun(@(f) f.hscale, F));
out.epslevel = max(cellfun(@(f) f.epslevel, F).*out.vscale)./max(out.vscale);

end
