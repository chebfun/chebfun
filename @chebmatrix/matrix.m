function varargout = matrix(A, varargin)
%MATRIX   Discretize a CHEBMATRIX as an ordinary matrix.
%   M = MATRIX(A, DIM) discretizes each block in the chebmatrix A using the
%   dimension vector DIM for all functions. In case the domain of A has
%   breakpoints, the vector DIM must specify the desired discretization
%   dimension for each subinterval.
%
%   MATRIX(A, DIM, DOMAIN) replaces the 'native' domain of A with DOMAIN.
%   Usually this would be done to introduce a breakpoint.
%
%   MATRIX(A,...,DISCTYPE) uses the chebDiscretization whose consructor is
%   DISCTYPE. The default is set by CHEBOPPREF. 
%
%   Example:
%     d = [0 1];
%     A = [ operatorBlock.eye(d), operatorBlock.diff(d) ];
%     matrix(A, 5, @chebcolloc2)
%     matrix(A, 5, @ultraS)
%
% See also CHEBOPPREF, CHEBDISCRETIZATION, CHEBDISCRETIZATION/MATRIX. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Any non-numeric argument should be a chebDiscretization constructor. 
discType = [];
numericargs = cellfun(@isnumeric, varargin);
for k = find( ~numericargs(:)' )
    discType = varargin{k};
    varargin(k) = [];  % Delete from list.
end

% Apply the default if needed.
if ( isempty(discType) || ~isa(discType, 'function_handle') )
    p = cheboppref;
    discType = p.discretization;
end

% Discretize.
d = discType(A, varargin{:});
[varargout{1:nargout}] = matrix(d);

end
