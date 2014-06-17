function out = matrix(N, varargin)
%   OUT = MATRIX(N, DIM) returns an DIM-point discretization of the linear
%   operator N. If N is not linear an error is thrown. 
%
%   OUT = MATRIX(N, DIM, 'oldschool') forces the returned differentiation
%   matrices to be square, rather than rectangular. See LINOP/MATRIX and
%   LINOP/FEVAL for further details.
%
% See also FEVAL, LINOP/MATRIX, LINOP/FEVAL.

[L, ignored, fail] = linop(N);
if ( fail )
    error('CHEBFUN:CHEBOP:matrix:nonlinear',...
        'Matrix expansion is only allowed for linear CHEBOP objects.')
end

if ( (numel(varargin) > 1) && strcmpi(varargin{2}, 'oldschool') )
    warnState = warning('off', 'CHEBFUN:LINOP:feval:deprecated');
    out = feval(L, varargin{:});
    warning(warnState);
else
    out = matrix(L, varargin{:});
end

end
