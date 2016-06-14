function f = constructor(f, op, varargin)
%CONSTRUCTOR   Main CHEBFUN3T constructor.
%   Classical full tensor approach for 3D functions. No low-rank technique
%   is involved here. This is experimental code, not for ordinary use.
%
% See also CHEBFUN3.  

%   The 'trig' flag is NOT implemented.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse the inputs:
[op, dom, pref, vectorize] = parseInputs(op, varargin{:});

% Set preferences:
tech        = pref.tech();
tpref       = tech.techPref;
grid = tpref.minSamples;
maxSample = 363;
% Since 363*363*363=(6916)^2, this is equivalent to working with a matrix 
% of size 6916 * 6916 which needs 380MB of memory.
pseudoLevel = pref.cheb3Prefs.chebfun3eps; 
prefStruct = pref.cheb2Prefs;
passSampleTest = prefStruct.sampleTest;

[xx, yy, zz] = ndgrid(dom(1:2), dom(3:4), dom(5:6));
A = op(xx, yy, zz);
if ( isscalar(A) )
    op = @(x,y,z) op(x,y,z) + 0*x + 0*y + 0*z;
end

%%
isHappy = 0;
m = grid; n = grid; p = grid;
for i=1:10
    %% Compute values at 3D chebpts
    out = tech.tensorGrid([m, n, p], dom);
    xx = out{1};
    yy = out{2};
    zz = out{3};
    F = op(xx,yy,zz);
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(F(:)));
    if ( isinf(vscale) )
        error('CHEBFUN:CHEBFUN3T:constructor:inf', ...
            'Function returned INF when evaluated');
    elseif ( any(isnan(F(:)) ) )
        error('CHEBFUN:CHEBFUN3T:constructor:nan', ...
            'Function returned NaN when evaluated');
    end
    
    %% Call 3D vals2coeffs
    coeffs3D = chebfun3t.vals2coeffs(F);
    % happinessCheck
    [isHappyX, isHappyY, isHappyZ, cutoffX2, cutoffY2, cutoffZ2] = ...
        happinessCheck3D(coeffs3D, pref);
    isHappy = isHappyX & isHappyY & isHappyZ;
    
    %% PHASE II
    % This phase refines the grid ONLY in "unhappy" directions.
    % We don't refine here if ALL the directions are nonhappy (as mentioned 
    % in the following). If, however, less than 3 directions are nonhappy, 
    % we refine in those directions only.
    
    m2 = m; n2 = n; p2 = p;
    while (~isHappy)
        % We need to refine the grid. We refine in 2 different ways: If all
        % 3 directions are unhappy, the new grid is bigger in each
        % direction by a factor of sqrt(2). If 1 or 2 of the directions are
        % unhappy, we refine (only those unhappy directions) by a factor of
        % 2.
        
        if (~isHappyX && ~isHappyY && ~isHappyZ)
            % All directions are unhappy. To avoid getting an unncecessarily
            % big tensor, refine differently.
            break;
        end
        if (~isHappyX)
            m2 = gridRefinePhase2(m2, pref);
        end
        if (~isHappyY)
            n2 = gridRefinePhase2(n2, pref);
        end
        if (~isHappyZ)
            p2 = gridRefinePhase2(p2, pref);
        end
        out = tech.tensorGrid([m2, n2, p2], dom);
        xx = out{1};
        yy = out{2};
        zz = out{3};
        F = op(xx,yy,zz);
        coeffs3D = chebfun3t.vals2coeffs(F);
        [isHappyX, isHappyY, isHappyZ, cutoffX2, cutoffY2, cutoffZ2] = ...
            happinessCheck3D(coeffs3D,pref);
        isHappy = isHappyX & isHappyY & isHappyZ;
    end
    
     if ( isHappy )
         % Chop the tail:
         coeffs3D = coeffs3D(1:cutoffX2, 1:cutoffY2, 1:cutoffZ2);
         
         % Construct a CHEBFUN3T object:
         f.coeffs = coeffs3D;
         f.vscale = max(abs(F(:)));
         f.domain = dom;
         
         % Step 2: Apply a SAMPLETEST:
         [m2, n2, p2] = size(coeffs3D);
         if ( passSampleTest )
             tol2 = getTol3D(xx, yy, zz, F, length(xx), dom, pseudoLevel);
             pass = sampleTest(f, op, tol2, vectorize);
             if ( pass ) % :)
                 break
             end 
         end
     end
    m = gridRefinePhase1(m);
    n = gridRefinePhase1(n);
    p = gridRefinePhase1(p);
end

% Reached maximum allowed tensor size and still unhappy?
if ( i == 10 && ~isHappy )
    % Throw a warning and construct from the latest unhappy approximation
    warning('CHEBFUN:CHEBFUN3T:constructor:notResolved', ...
                'Unresolved with maximum CHEBFUN3T length: %u.', maxSample);
    f.coeffs = coeffs3D;
    f.vscale = max(abs(F(:)));
    f.domain = dom;
end

end

%%
function [isHappyX, isHappyY, isHappyZ, cutoffX2, cutoffY2, cutoffZ2] = ...
    happinessCheck3D(simple_3D_coeffs, pref)
% An alternative technique: check whether all the rows, columns and tubes have been resolved.
% Somewhat similar to Chebfun2 constructor (but with all the rows, cols and tubes instead of the pivot ones)

% TODO: THIS DOES NOT PROPERLY WORK WITH TRIGs.

colChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 1), 2);
fCol = chebtech2({[], colChebtech});
[isHappyX, cutoffX2] = happinessCheck(fCol, [], [], [], pref);
    
rowChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 2), 2);
fRow = chebtech2({[], rowChebtech});
[isHappyY, cutoffY2] = happinessCheck(fRow, [], [], [], pref);

tubeChebtech = sum(chebfun3t.unfold(abs(simple_3D_coeffs), 3), 2);
fTube = chebtech2({[], tubeChebtech});
[isHappyZ, cutoffZ2] = happinessCheck(fTube, [], [], [], pref);
end
%%
function grid = gridRefinePhase1(grid)
% Grid refinement strategy in Phase 1. It does the same for all TECHS.
grid = floor(sqrt(2)^(floor(2*log2(grid)) + 1)) + 1;
end

%%
function grid = gridRefinePhase2(grid, pref)
% Grid refinement strategy for tech in Phase 2 

% What tech am I based on?:
tech = pref.tech();

% What is the next grid size?
if ( isa(tech, 'chebtech2') )
    % Double sampling on tensor grid:
    grid = 2*grid-1;
elseif ( isa(tech, 'trigtech') )
    % Double sampling on tensor grid:
    grid = 2^(floor(log2(grid)+1));
elseif ( isa(tech, 'chebtech1') )
    grid = 3 * grid;
else
    error('CHEBFUN:CHEBFUN3:constructor:gridRefinePhase2:techType', ...
        'Technology is unrecognized.');
end

end
%%

function tol = getTol3D(xx, yy, zz, vals, grid, dom, pseudoLevel)
%See https://github.com/chebfun/chebfun/issues/1491
[m,n,p] = size(vals);
% Remove some edge values so that df_dx, df_dy and df_dz have the same size. 
df_dx = diff(vals(:,1:n-1,1:p-1),1,1) ./ diff(xx(:,1:n-1,1:p-1),1,1); % xx changes in the first mode.
df_dy = diff(vals(1:m-1,:,1:p-1),1,2) ./ diff(yy(1:m-1,:,1:p-1),1,2); % yy changes row-wise (2nd mode).
df_dz = diff(vals(1:m-1,1:n-1,:),1,3) ./ diff(zz(1:m-1,1:n-1,:),1,3); % zz changes tube-wise (3rd mode).

J = max(max(abs(df_dx),abs(df_dy)), abs(df_dz));
Jac_norm = max(J(:)); % An approximation for the norm of the Jacobian over the whole domain.
vscale = max(abs(vals(:)));
tol = grid^(4/5) * max( abs(dom(:) ) ) * max(Jac_norm, vscale) * pseudoLevel;
end
%%
function [op, dom, pref, vectorize] = parseInputs(op, varargin)
vectorize = 0;
pref = chebfunpref();

if ( isa(op, 'char') )     % CHEBFUN3T( CHAR )
    op = str2op(op);
end

if nargin == 1
    dom = [-1 1 -1 1 -1 1];
elseif ( ( (nargin == 2) && ~any(strcmpi(varargin{1}, {'trig', 'periodic'})) ))    
    % DOMAIN is specified
    dom = varargin{1};
elseif ( ( (nargin == 2) && any(strcmpi(varargin{1}, {'trig', 'periodic'})) ))    
    % TECH is specified
    pref.tech = @trigtech;    
    dom = [-1 1 -1 1 -1 1];
elseif ( ( (nargin == 3) && any(strcmpi(varargin{2}, {'trig', 'periodic'})) ))
    % DOMAIN and TECH are specified
    dom = varargin{1};
    pref.tech = @trigtech;
elseif ( ( (nargin == 3) && any(strcmpi(varargin{1}, {'trig', 'periodic'})) ))
    % TECH and DOMAIN are specified
    pref.tech = @trigtech;
    dom = varargin{3};
elseif ( (nargin == 3) && strcmpi(varargin{1}, 'eps') )
    % EPS is specified
    pref.chebfuneps = varargin{2};
    dom = [-1 1 -1 1 -1 1];
elseif ( (nargin == 4) && strcmpi(varargin{1}, 'eps') ...
        && any(strcmpi(varargin{3}, {'trig', 'periodic'})) )
    % EPS and TECH are specified
    pref.tech = @trigtech;
    pref.chebfuneps = varargin{3};
    dom = [-1 1 -1 1 -1 1];
elseif ( (nargin == 4) && strcmpi(varargin{2}, 'eps') && ...
        any(strcmpi(varargin{1}, {'trig', 'periodic'})) )
    % TECH and EPS are specified
    pref.tech = @trigtech;
    pref.chebfuneps = varargin{3};
    dom = [-1 1 -1 1 -1 1];
elseif ( (nargin == 4) && strcmpi(varargin{2}, 'eps') && ...
        ~any(strcmpi(varargin{1}, {'trig', 'periodic'})) )
    % DOMAIN and EPS are specified
    pref.chebfuneps = varargin{3};
    dom = varargin{1};    
elseif ( (nargin == 4) && strcmpi(varargin{1}, 'eps') && ...
        ~any(strcmpi(varargin{3}, {'trig', 'periodic'})) )
    % EPS and DOMAIN are specified
    pref.chebfuneps = varargin{2};
    dom = varargin{3};
elseif ( (nargin == 5) && any(strcmpi(varargin{2}, {'trig', 'periodic'}))...
        &&  strcmpi(varargin{3}, 'eps') )
    % DOMAIN, TECH and EPS are specified
    dom = varargin{1};
    pref.tech = @trigtech;
    pref.chebfuneps = varargin{4};
elseif ( (nargin == 5) && strcmpi(varargin{2}, 'eps') ...
        && any(strcmpi(varargin{4}, {'trig', 'periodic'})) )
    % DOMAIN, EPS and TECH are specified
    dom = varargin{1};
    pref.chebfuneps = varargin{3};
    pref.tech = @trigtech;
elseif ( (nargin == 5) && any(strcmpi(varargin{1}, {'trig', 'periodic'}))...
        &&  strcmpi(varargin{2}, 'eps') )
    % TECH, EPS and DOMAIN are specified
    pref.tech = @trigtech;    
    pref.chebfuneps = varargin{3};
    dom = varargin{4};
elseif ( (nargin == 5) && any(strcmpi(varargin{1}, {'trig', 'periodic'}))...
        &&  strcmpi(varargin{3}, 'eps') )
    % TECH, DOMAIN and EPS are specified
    pref.tech = @trigtech;    
    dom = varargin{2}; 
    pref.chebfuneps = varargin{4};
elseif ( (nargin == 5) && any(strcmpi(varargin{4}, {'trig', 'periodic'}))...
        &&  strcmpi(varargin{1}, 'eps') )
    % EPS, DOMAIN and TECH are specified
    pref.chebfuneps = varargin{2};
    dom = varargin{3};
    pref.tech = @trigtech;    
elseif ( (nargin == 5) && any(strcmpi(varargin{3}, {'trig', 'periodic'}))...
        &&  strcmpi(varargin{1}, 'eps') )
    % EPS, TECH and DOMAIN are specified
    pref.chebfuneps = varargin{2};
    pref.tech = @trigtech;
    dom = varargin{4};
else
    error('CHEBFUN:CHEBFUN3T:constructor', ...
                    'Domain not valid or fully determined.');
end
end
%%
function op = str2op( op )
% OP = STR2OP(OP), finds independent variables in a string and returns an op
% handle than can be evaluated.

vars = symvar(op); % Independent variables
if ( numel(vars) > 3)
    error('CHEBFUN:CHEBFUN3T:constructor:str2op:depvars', ...
        'Too many independent variables in string input.');
else
    op = eval(['@(' vars{1} ',' vars{2} ',' vars{3} ')' op]);
end

end