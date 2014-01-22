function M = matrix(A,varargin)
%MATRIX  Discretize a chebmatrix as an ordinary matrix.
%   M = MATRIX(A,DIM) discretizes each block in the chebmatrix A using the
%   dimension vector DIM for all functions.
%
%   M = MATRIX(A,DIM,DOMAIN) replaces the 'native' domain of A with DOMAIN.
%   Usually this would be done to introduce a breakpoint.
%
%   The discretization type is controlled by the 'prefs' property of A.
%
%   Example:
%     d = [0 1];
%     A = [ operatorBlock.eye(d), operatorBlock.diff(d) ];
%     A.prefs.discretization = @colloc2;
%     matrix(A,5)
%     A.prefs.discretization = @ultraS;
%     matrix(A,5)
%
%   See also CHEBOPPREF, CHEBDISCRETIZATION. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

d = A.prefs.discretization(A,varargin{:});
M = matrix(d);

end