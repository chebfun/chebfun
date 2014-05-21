function G = quasimatrix(F, varargin)
%QUASIMATRIX   A quasimatrix is an array of CHEBFUN objects.
%
%   TODO: Use this to document the difference between an array-valued CHEBFUN
%   and a quasimatrix.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% This method is simply a wrapper for CHEB2QUASI() of the CHEBFUN constructor:

% TODO: Document this file.

if ( iscell(F) && any(cellfun(@(f) isa(f, 'chebfun'), F)) )
    % Construct from a cell of CHEBFUNs.
    G = chebfun.cell2quasi(F);
    return
end

if ( ~isa(F, 'chebfun') )
    % Create a CHEBFUN:
    F = chebfun(F, varargin{:});
end

% Convert to a quasimatrix:
G = cheb2quasi(F);

end
