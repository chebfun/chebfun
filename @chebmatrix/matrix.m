function M = matrix(A, varargin)
%MATRIX  Discretize a chebmatrix as an ordinary matrix.
%   M = MATRIX(A, DIM) discretizes each block in the chebmatrix A using the
%   dimension vector DIM for all functions. In case the domain of A has
%   breakpoints, the vector DIM must specify the desired discretisation
%   dimension for each subinterval.
%
%   M = MATRIX(A, DIM, DOMAIN) replaces the 'native' domain of A with DOMAIN.
%   Usually this would be done to introduce a breakpoint.
%
%   The discretization type is controlled by the 'prefs' property of A.
%
%   Example:
%     d = [0 1];
%     A = [ operatorBlock.eye(d), operatorBlock.diff(d) ];
%     A.prefs.discretization = @colloc2;
%     matrix(A, 5)
%     A.prefs.discretization = @ultraS;
%     matrix(A, 5)
%
%   See also CHEBOPPREF, CHEBDISCRETIZATION, CHEBDISCRETIZATION/MATRIX. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% TODO: Confirm with TAD why this was A.prefs.discretizatino, calling a
% dependent property? AB, 27/1/14.
d = A.prefs.discretisation(A, varargin{:});
M = matrix(d);

end