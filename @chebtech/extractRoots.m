function [g, n] = extractRoots(f, m, sides, tol)
% EXTRACTROOTS   Extract roots from ends of chebtechs.
%   EXTRACTROOTS(F) tries to extract roots of F with integer multiplicity 
%   at -1 and 1, so that the multiplicities of the remaining roots at -1 
%   and 1 are less than 1. The output G is a CHEBTECH and N specifies the 
%   multiplicity of the roots which have been extracted.   
%
%   EXTRACTROOT(F, M) extracts roots of F with multiplicity M at -1 and 1.
%   M and N should be identical.
%
%   EXTRACTROOT(F, M, SIDES) extracts roots of F with multiplicity M at the
%   endpoint(s) specified by SIDES, -1 or 1 or both. Note that the
%   dimension of M and that of SIDES should match.
%
%   EXTRACTROOT(F, [], SIDES) extracts roots of F at the specified
%   endpoints, so that the multiplicities of the remaining roots at the
%   endpoints are less than 1.
%
%   EXTRACTROOT(F, M, SIDES, TOL) is same as when TOL is omitted. But it 
%   uses TOL to determine the root multiplicity, so that G(-1) > TOL and
%   G(1) > TOL. Note that M or SIDES or both can be omitted.
%
% See also REPLACEROOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

if nargin < 4
    tol = 100*f.epslevel;
end

if ( nargin < 3 ) || ( isempty(sides) )
    sides = ones(1, 2); 
end

if ( nargin < 2 ) || ( isempty(m) )
    m = Inf(1, 2); 
end

valEnd = abs(feval(f, [-1 1]));
n = 0;

persistent D

while ( ( sides > 0 ) && ( valEnd < tol ) && ( all(n < m) ) )    
    
    coeffsF = flipud(f.coeffs);
    
    % Construct the matrix for the recurrence.
    n = length(f);
    v = ones(n-2,1);
    DL = spdiags([v 2*v v], 0:2, n-1, n-1);
    DL(1,1) = 1;
    
    % Construct G:
    g = f.make();
    
    % Computing the coefficients of G by solving the recurrence system.
    coeffsG = flipud(DL\coeffsF(2:end));
    
    % Compute the values at Chebyshev points.
    g.vals = f.chebpolyval(coeffsF);
    g.coeffs = coeffsG;
    g.vscale = max(g.vals);
    g.hscale = 1;
    g.ishappy = 1;
    g.epslevel = g.epslevel;  %[TODO]: This need to be fixed.
    
    % Update f0
    valEnd = abs(vals([1 end]));
    valEnd(~sides) = inf;
    num = num+1;

end

persistent DL
persistent DR

while ( n < m ) || all( valEnd(sides) )
    
    coeffsF = flipud(f.coeffs);
    
    % Construct the matrix for the recurrence.
    n = length(f);
    v = 0.5*ones(n-1,1);
    DL = spdiags([v 2*v v], 0:2, n-1, n-1);
    DL(1,1) = 1;
    
    % Construct G:
    g = f.make();
    
    % Computing the coefficients of G by solving the recurrence system.
    coeffsG = flipud(DL\coeffsF(2:end));
    
    % Compute the values at Chebyshev points.
    g.vals = f.chebpolyval(coeffsF);
    g.coeffs = coeffsG;
    g.vscale = max(g.vals);
    g.hscale = 1;
    g.ishappy = 1;
    g.epslevel = g.epslevel; %[TODO]: This need to be fixed.
    
    % Update f0
    valEnd = abs(vals([1 end]));
    valEnd(~sides) = inf;
    num = num+1;
    
end

