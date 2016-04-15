function f = constructor(f, op, varargin)
%CONSTRUCTOR   The main CHEBFUN3T constructor.
%   Classiacal full tensor approach for 3D functions. No low-rank technique
%   is involved here.

%TODO: add `trig` to chebfun3t: HAPPINESSCHECK here DOES NOT WORK WITH TRIGs PROPERLY.

% Parse the inputs:
[op, dom, pref, vectorize] = parseInputs(op, varargin{:});
% TODO: Make parseInputs the same as CHEBFUN3.

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
    [xx, yy, zz] = points3D(m, n, p, dom, pref);
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
    % This phase refines the grid ONLY in the directions which are needed. 
    % We don't refine here if ALL the directions are nonhappy (as mentioned 
    % in the following). If, however, 1 or 2 of the directions are nonhappy, 
    % we refine in those directions only.
    
    % TODO: Currently, in this Phase II, we allow the lengths of the grid 
    % in the directions which are not happy to be ONLY 2 times bigger than
    % the "happy directions". In other words, we do refinement only in ONE
    % LEVEL. Consider e.g., the following "unsymmetric" function:
    % f = @(x,y,z) sin(x./2 + y./3 + 20*z);
    % Here, f3t = chebfun3t(f) will be happy in X and Y directions in the 
    % first run of Phase I, i.e., with m = 12 and n = 12 but with
    % isHappyZ=0 when p = 12. What we do currently in this Phase II, is
    % that we refine in the Z direction only ONCE, i.e., we increase P to
    % 17 and check happiness3D again. If we are not still happy, then we go
    % to the next run of Phase I, i.e., we set m=n=p=35 and start again.
    % The point is that after we ran Phase II once, we check
    % happiness and if we were not happy in the same directions that we
    % were not BEFORE Phase II, we might allow refinement in those unhappy
    % directions ONLY AGAIN.
    % It is correct that the final chopping will reduce the length in the
    % directions we were happy before, the point is that with this present
    % constructor, we might increase the lengths even in the direction we
    % don't really need which "might" (depending on how big are m, n and p already) 
    % make the computations slow without any good reasons.
    
    m2 = m; n2 = n; p2 = p;
    while (~isHappy)
        % We need to refine the grid. We refine in 2 different ways: If all
        % 3 directions are unhappy, the new grid is bigger in each
        % direction by a factor of sqrt(2). If 1 or 2 of the directions are
        % unhappy, we refine only those unhappy directions by a factor of
        % 2.
        
        if (~isHappyX && ~isHappyY && ~isHappyZ)
            % All directions are unhappy. Refine differently to avoid
            % rapidly getting to a huge tensor.
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
        [xx, yy, zz] = points3D(m2, n2, p2, dom, pref);
        F = op(xx,yy,zz);
        coeffs3D = chebfun3t.vals2coeffs(F);
        [isHappyX, isHappyY, isHappyZ, cutoffX2, cutoffY2, cutoffZ2] = ...
            happinessCheck3D(coeffs3D,pref);
        isHappy = isHappyX & isHappyY & isHappyZ;
    end
    
     if ( isHappy )
         % Chop the tail.
         coeffs3D = coeffs3D(1:cutoffX2, 1:cutoffY2, 1:cutoffZ2);
         
         % Construct a CHEBFUN3T object:
         f.coeffs = coeffs3D;
         f.vscale = max(abs(F(:)));
         f.domain = dom;
         
         % Step 2: Start a SAMPLETEST
         [m2, n2, p2] = size(coeffs3D);
         
         if ( passSampleTest )
             tol2 = getTol3D(xx, yy, zz, F, length(xx), dom, pseudoLevel);
             pass = sampleTest(f, op, tol2, vectorize);
             if ( pass ) % Step 1 was happy and sampleTest is also passed.
                 break
             end 
         end
     end
    m = gridRefinePhase1(m);
    n = gridRefinePhase1(n);
    p = gridRefinePhase1(p);
end

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
function [xx, yy, zz] = points3D(m, n, p, dom, pref)
% Get the sample points that correspond to the right grid for a particular
% technology.

% What tech am I based on?:
tech = pref.tech();

if ( isa(tech, 'chebtech2') )
    x = chebpts( m, dom(1:2), 2 );   % x grid.
    y = chebpts( n, dom(3:4), 2 );
    z = chebpts( p, dom(5:6), 2 );
    [xx, yy, zz] = ndgrid( x, y, z ); 
elseif ( isa(tech, 'chebtech1') )
    x = chebpts( m, dom(1:2), 1 );   % x grid.
    y = chebpts( n, dom(3:4), 1 ); 
    z = chebpts( p, dom(5:6), 1 );
    [xx, yy, zz] = ndgrid( x, y, z ); 
elseif ( isa(tech, 'trigtech') )
    x = trigpts( m, dom(1:2) );   % x grid.
    y = trigpts( n, dom(3:4) );
    z = chebpts( p, dom(5:6), 2 );
    [xx, yy, zz] = ndgrid( x, y, z );
else
    error('CHEBFUN:CHEBFUN3T:constructor:points3D:tecType', ...
        'Unrecognized technology');
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
    grid = 2^(floor( log2(grid)+1));
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