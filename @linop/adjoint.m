function [Lstar, op, bcOpL, bcOpR, bcOpM] = adjoint(L, bcType)
%ADJOINT   Compute the adjoint of a LINOP.
%   ADJOINT(L), where L is a LINOP, returns the adjoint LINOP of L under
%   the assumption that L only has endpoint or periodic functional constraints.
%
%   ADJOINT(L, BCTYPE) allows for more general boundary conditions.
%
% See also ?.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

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
% Create adjoint linop
% periodic boundary conditions
if ( strcmp(bcType, 'periodic') )
    bcOpL = [];
    bcOpR = [];
    bcOpM = [];
    Lstar.constraint = L.constraint;
% general boundary conditions
else
    % Adjoint boundary conditions
    [Cstar, bcOpL, bcOpR, bcOpM] = adjointBCs(L, bcType);
    Lstar.constraint = Cstar;
end

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
