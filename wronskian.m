function w = wronskian(L, varargin)
%WRONSKIAN   Wronskian of chebfuns.
%   WRONSKIAN(L, f1, ..., fn) computes the wronskian of the chebfuns.
%   L is the linear differential operator and f1, ..., fn are solutions of the
%   homogenous problem. The algorithm is based on Abel's identity for computing
%   the wronskian.
%   
%   WRONSKIAN(L, F) does the same where F is a quasimatrix or an array valued
%   chebfun or a chebmatrix. If F is a chebmatrix, then the form of F is assumed
%   to be a chebmatrix with a single column with n chebfuns.
%
% See also LINOP

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = {};
% Quasimatrix case:
if ( isa(varargin{1}, 'cell') )
    F = varargin{1};
else
    if ( isa(varargin{1}, 'chebfun') )
        G = varargin{1};
        nCols = size(G, 2);
        % Array valued chebfuns
        if (  nCols > 1 )
            for i = 1:nCols
                F{i} = G(:, i); %#ok<AGROW>
            end
        else
            % Chebfuns:
            for i = 1:length(varargin)
                F{i} = varargin{i}; %#ok<AGROW>
            end
        end
    else
        if ( isa(varargin{1}, 'chebmatrix') )
            G = varargin{1};
            G = (G.blocks).';
            w = wronskian(L, G);
            return                
        end
    end
end
% Extract the coefficients of the operator:
L = linop(L);
c = toCoeff(L.blocks{1});
% The (n-1)th coefficient in standard form:
p = c{2}./c{1};
dom = L.domain;
a = dom(1);
n = size(c, 2)-1;
nFuns = size(F, 2);
if ( nFuns ~= n )
    error( 'CHEBFUN:WRONSKIAN', 'Number of chebfuns is not the same as the order of the operator' )
end

% [TODO]: This can be vectorized using fancy chembatrices etc?
W = zeros(n);
for i = 1:n
    for j = 1:n
        f = F{j};
        W(i, j) = f(a);
        F{j} = diff(f); %#ok<AGROW>
    end   
end

% Compute the determinant at the left end of the domain:
A = det(W);

% Apply Abel's identity
w = A*exp(-cumsum(p));
% Make sure w(a) = A
w = (w - w(a)) + A;
end
    
    