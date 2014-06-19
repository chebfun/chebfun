function varargout = polyeigs(varargin)
%POLYEIGS   Polynomial CHEBOP eigenvalue problem.
%   POLYEIGS of a CHEBOP is not yet supported for p > 1.
%
% See also EIGS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Find CHEBOP inputs:
isChebop = cellfun(@(v) isa(v, 'chebop'));

if ( nargin > 2 && any(isChebop(3:end)) )
    
    % POLYEIGS is not yet supported.
    error('CHEBFUN:CHEBOP:polyeigs:notSupported', ...
        'CHEBOP/POLYEIGS() is not currently supported.');
    
else
    
    % Simply call EIGS for this problem:
    [varargout{1:nargout}] = eigs(varargin{:});
    
end

end
