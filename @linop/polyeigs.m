function varargout = polyeigs(varargin)
%POLYEIGS   Polynomial LINOP eigenvalue problem.
%   POLYEIGS of a LINOP is not yet supported for p > 1.
%
% See also EIGS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Find LINOP inputs:
isLinop = cellfun(@(v) isa(v, 'linop'));

if ( nargin > 2 && any(isLinop(3:end)) )
    
    % POLYEIGS is not yet supported.
    error('CHEBFUN:LINOP:polyeigs:noSupport', ...
        'CHEBOP/POLYEIGS() is not currently supported.');
    
else
    
    % Simply call EIGS for this problem:
    [varargout{1:nargout}] = eigs(varargin{:});
    
end

end
