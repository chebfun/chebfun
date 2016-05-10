function eig(varargin)
%   EIG(N) is not supported. Use EIGS(N) to find selected eigenvalues of a
%   linear operator.
%
% See also CHEBOP/EIGS

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBOP:eig:useEIGS',...
    'Use EIGS to find selected eigenvalues of a linear operator.')

end
