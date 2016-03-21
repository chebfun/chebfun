function [Lstar, op, bcOpL, bcOpR, bcOpM] = adjoint(L, bcType)
%ADJOINT   Compute the adjoint of a LINOP.
%   [LSTAR, OP, BCOPL, BCOPR, BCOPM] = ADJOINT(L, BCTYPE) computes the adjoint
%   of the LINOP L. ADJOINT requires that L represents a linear differential
%   operator with either endpoint of periodic boundary conditions. If L
%   represents a system of differential equations then the highest order
%   derivative in each variable must be the same. Integral operators and exotic
%   functional constraints are not supported.
%
%   The output is a LINOP LSTAR that represents the adjoint differential
%   operator and adjoint constraints, as well as four function handles that can
%   be used to construct a CHEBOP. OP is a function handle for the differential
%   operator and BCOPL, BCOPR and BCOPM are function handles for the left,
%   right and mixed boundary conditions respectively. 
%
% See also ADJOINTFORMAL and ADJOINTBCS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% 
% Check for 2 inputs
if ( nargin < 2 )
    error('CHEBFUN:LINOP:adjoint:inputs', ...
    'ADJOINT requires two inputs.')
end

%% 
% Check to so if we can compute an adjoint for this operator!
parseInputs(L, bcType);

%%
% Formal adjoint:
pref = chebfunpref();
if ( strcmpi(bcType, 'periodic') )
    pref.tech = @trigtech;
end 
[Lstar, op] = adjointFormal(L, pref);

% Trivial case:
if ( max(max(L.diffOrder)) == 0 )
    bcOpL = [];
    bcOpR = [];
    bcOpM = [];
    return
end

%%
% adjoint boundary conditions
[Cstar, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType);
Lstar.constraint = Cstar;

end



function parseInputs(L, bcType)
% function to parse inputs and catch errors

[m,n] = size(L);

if ( min(min(L.diffOrder)) < 0 )
    % Check for integral operators:
    error('CHEBFUN:LINOP:adjoint:difforder', ...
    'ADJOINT doesn''t support integral operators for the moment.')
elseif ( size(L.domain, 2) > 2 )
    % [TODO]: Support piecewise domains.
    error('CHEBFUN:LINOP:adjoint:domain', ...
        'ADJOINT doesn''t support piecewise domains for the moment.');
elseif ( ~any(strcmp(bcType, {'periodic', 'bvp', 'ivp', 'fvp'})) )
    error('CHEBFUN:LINOP:adjoint:boundaryconditions', ...
    'ADJOINT doesn''t support this type of boundary conditions for the moment.');
elseif ( m ~= n )
    % [TODO]: Support nonsquare systems.
    error('CHEBFUN:LINOP:adjoint:systems', ...
    ['ADJOINT doesn''t support nonsquare systems at the moment.']);
elseif ( n > 1 && any(diff(max(L.diffOrder)) ~= 0) )
    % [TODO]: Support block systems with arbitrary differential order in 
    %         each variable.
    error('CHEBFUN:LINOP:adjoint:systems', ...
    ['ADJOINT doesn''t support systems with arbitrary ',...
      'differential order in each variable at the moment.']);
end

end
