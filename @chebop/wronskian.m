function w = wronskian(N, varargin)
%WRONSKIAN   Wronskian of chebfuns.
%   WRONSKIAN(N, f1, ..., fn) computes the wronskian of the CHEBFUN objects fk.
%   N is the linear differential operator and f1, ..., fn are solutions of the
%   homogenous problem. The algorithm is based on Abel's identity for computing
%   the wronskian.
%   
%   WRONSKIAN(N, F) does the same where F is a quasimatrix or an array-valued
%   CHEBFUN or a CHEBMATRIX. If F is a CHEBMATRIX, then the form of F is assumed
%   to be a CHEBMATRIX with a single column with n CHEBFUNs.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert to a Linop (if possible):
[L, ignored, fail] = linop(N); %#ok<ASGLU>
if ( fail ) % Throw an error if the operator is not linear.
    error('CHEBFUN:CHEBOP:wronskian:nonlinear', ...
        'WRONSKIAN() only support linear operators.')
end

% Call LINOP/WRONSKIAN():
w = wronskian(L, varargin{:});

end