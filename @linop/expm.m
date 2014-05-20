function u = expm(L, t, u0, prefs)
%EXPM      Exponential semigroup of a LINOP.
%   u = EXPM(L, T, U0) uses matrix exponentiation to propagate an initial
%   condition U0 for time T through the differential equation u' = L*u, where L
%   is a LINOP. Formally, the solution is given by u(t) = exp(t*L)*u0, where
%   exp(t*L) is a semigroup generated by L.
%
%   If T is a vector, then U is a CHEBMATRIX and will have one column for each
%   entry of T. If T is a scalar and U0 is a scalar CHEBFUN or CHEBMATRIX, then
%   U is a CHEBFUN. Note that this latter behaviour differs from other LINOP
%   methods, which usually return a CHEBMATRIX even for the scalar case.
%
%   L should have appropriate boundary conditions to make the problem
%   well-posed. Those conditions have zero values; i.e. are represented by
%   B*u(t)=0 for a linear functional B.
%
%   EXPM(..., PREFS) accepts a preference structure or object like that created
%   by CHEBOPPREF.
%
%   EXAMPLE: Heat equation
%      d = [-1 1];  x = chebfun('x', d);
%      D = operatorBlock.diff(d);  
%      A = linop( D^2 );  
%      E = functionalBlock.eval(d);
%      A = addBC(A, E(d(1)), 0);   % left Dirichlet condition
%      A = addBC(A, E(d(2)), 0);   % right Dirichlet condition
%      u0 = exp(-20*(x+0.3).^2);  
%      t = [0 0.001 0.01 0.1 0.5 1];
%      u = expm(A, t, u0);
%      colr = zeros(6, 3);  colr(:,1) = 0.85.^(0:5)';
%      clf, set(gcf, 'defaultaxescolororder', colr)
%      plot(chebfun(u), 'linewidth', 2)
%
% See also LINOP/ADDBC.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 4 )
    prefs = cheboppref;
end

discType = prefs.discretization;
isFun = isFunVariable(L); 

%% Set up the discretization:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    disc = discType(L);  
    
    % Merge domains of the operator and the initial condition.
    disc.domain = chebfun.mergeDomains(disc.domain, u0.domain); 
    
    % Set the allowed discretisation lengths: 
    dimVals = prefs.dimensionValues;
    
    dimVals( dimVals < length(u0) ) = [];
    
    % Apply the discretiztion dimension on all pieces:
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
else
    % A discretization is given:
    disc = discType;
        
    % Initialise dimVals;
    dimVals = max(disc.dimension);
end

if ( isempty(L.continuity) )
     % Apply continuity conditions:
     disc.source = deriveContinuity(disc.source,disc.domain);
end

% Initialise happiness:
numInt = disc.numIntervals;
isDone = false(1, numInt);

if ( isa(u0, 'chebfun') )
    u0 = chebmatrix({u0}); 
elseif ( ~isa(u0, 'chebmatrix') )
    error('CHEBFUN:linop:expm:unknown', ...
        'No support for inputs of type %s.', class(u0));
end

%% Loop over different times.
allu = chebmatrix({});
for i = 1:length(t)
    
    %% Loop over a finer and finer grid until happy:
    for dim = dimVals
 
        disc.dimension(~isDone) = dim;

        % Discretize the operator (incl. constraints/continuity):
        [E, P] = expm(disc, t(i));
        
        % Discretize the initial condition.
        discu = disc;
        do = max(getDiffOrder(disc.source), 0);
        do = max(do, [], 1);
        for k = 1:numel(u0.blocks)
            discu.dimension = disc.dimension + +do(k);
            xIn = functionPoints(discu);
            if ( ~isnumeric(u0.blocks{k}) )
                f.blocks{k} = feval(u0.blocks{k}, xIn);
            end
        end
        v0 = cell2mat(f.blocks);  
        
        % Propagate.
        v = P*(E*v0);
        
        % Convert the different components into cells
        u = partition(disc, v);
        uFun = u(isFun);
        scale = max( cellfun(@max,uFun) );
        
        % Test the happieness of the function pieces:
        [isDone, epsLevel] = testConvergence(disc, uFun, scale, prefs);
        
        if ( all(isDone) )
            break
        end
        
    end
    
    if ( ~all(isDone) )
        warning('LINOP:expm:NoConverge', ...
            'Matrix exponential may not have converged.')
    end
    
    %% Tidy the solution for output:
    ucell = mat2fun(disc, u);
    doSimplify = @(f) simplify( f, max(eps, epsLevel) );
    ucell = cellfun( doSimplify, ucell, 'uniform', false );
    allu = [ allu, chebmatrix(ucell) ];
end

u = allu;

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( all(size(u) == [1 1]) )
    u = u{1};
end

end
