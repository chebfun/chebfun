function f = constructor( f, op, varargin )
%CONSTRUCTOR   BALLFUN constructor.
%   Given a function OP of three variables, defined on the unit ball, this
%   code represents it as a BALLFUN object. A BALLFUN object is a
%   tensor-product representation of a function that has been "doubled-up"
%   in the r- and theta-variables. This doubled-up function is represented
%   using in a Chebyshev-Fourier-Fourier basis.
%
% See also SPHEREFUN, DISKFUN, CHEBFUN2, CHEBFUN3.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[op, pref, isVectorize] = parseInputs(op, varargin{:});

% Set preferences:
tech            = pref.tech();
tpref           = tech.techPref;
grid1           = tpref.minSamples;
grid2           = tpref.minSamples;
grid3           = tpref.minSamples;
maxSample       = tpref.maxLength; % maxSample = max grid dimensions.

if ( isa(op, 'ballfun') )     % BALLFUN( BALLFUN )
    f = op;
    return
elseif ( isa(op, 'double') )   % BALLFUN( DOUBLE )
    f = constructFromDouble(f, op);
    return
end

isHappy = 0;     % We are currently unresolved.
failure = 0;     % Reached max discretization size without being happy.

while ( ~isHappy && ~failure )
    %% Main loop of the constructor
    vals = evaluate(op, [grid1, grid2, grid3], vectorize);
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('CHEBFUN:CHEBFUN3:constructor:inf', ...
            'Function returned INF when evaluated');
    elseif ( any(isnan(vals(:))) )
        error('CHEBFUN:CHEBFUN3:constructor:nan', ...
            'Function returned NaN when evaluated');
    end
    
    % If the rank of the function is above maxRank then stop.
    if ( max([grid1*grid2, grid2*grid3, grid1*grid3]) > maxSample )
        warning('CHEBFUN:BALLFUN:constructor:dimensions', ...
            'Not well-approximated by a Chebyshev-Fourier-Fourier expansion.');
        failure = 1;
    end
    
    [grid1, grid2, grid3, isHappy] = ballfunHappiness( vals, pref );
    
end

% We are now happy so make a BALLFUN from its values: 
f.coeffs = ballfun.vals2coeffs(vals);
end

%%
function f = constructFromDouble( f, op )
%CONSTRUCTFROMDOUBLE  Constructor BALLFUN from matrix of values.

if ( numel(op) == 1 )
    f = ballfun(@(x,y,z) op + 0*x, dom);
    return
end
f.coeffs = ballfun.vals2coeffs(op);
end


function vals = evaluate(g, S, isVectorize)
%EVALUATE   Evaluate at a Cheb-Fourier-Fourier grid of size S.
%  EVALUATE(g, S) returns the S(1)xS(2)xS(3) values of g at a
%  Chebyshev-Fourier-Fourier grid for the function
%  g(r, lambda, theta).

% Convert a handle_function to a ballfun function
m = S(1);
n = S(2);
p = S(3);

% Build the grid of evaluation points
r = chebpts(m);
lam = pi*trigpts(n);
th = pi*trigpts(p);

[rr, ll, tt] = ndgrid(r, lam, th);
if ~isVectorize
    % Evaluate function handle at tensor grid:
    vals = feval(g, rr, ll, tt);
else
    % If vectorize flag is turned out, then FOR loop: 
    vals = zeros(size(rr)); 
    for i1 = 1:size(vals,1)
        for j1 = 1:size(vals,2)
            for k1 = 1:size(vals,3)
                vals(i1,j1,k1) = g( rr(i1,j1,k1), ll(i1,j1,k1), tt(i1,j1,k1) ); 
            end
        end
    end
end

% Test if the function is constant
if size(vals) == 1
    vals = vals(1)*ones(S);
end

end

%%
function [op, pref, isVectorize] = parseInputs(op, varargin)
% Parse user inputs to BALLFUN.

isVectorize = 0;
isCoeffs = 0;
pref = chebfunpref();

% Preferences structure given?
isPref = find(cellfun(@(p) isa(p, 'chebfunpref'), varargin));
if ( any(isPref) )
    pref = varargin{isPref};
    varargin(isPref) = [];
end

if ( isa(op, 'char') )     % CHEBFUN3( CHAR )
    op = str2op(op);
end

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    elseif ( isnumeric(varargin{k}) )
        if ( numel(varargin{k}) == 3 ) % length is specified.
            % Interpret this as the user wants a fixed degree ballfun.
            S = varargin{k};
            op = evaluate(op, S, isVectorize);
        end
    elseif any(strcmpi(varargin{k}, {'vectorize', 'vectorise'}))
        isVectorize = true;
    elseif strcmpi(varargin{k}, 'coeffs')
        isCoeffs = 1;
    end
end

if ( isCoeffs )
    op = ballfun.coeffs2vals(op);
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~isVectorize && ~isnumeric(op) ) % another check
    [isVectorize, op] = vectorCheck(op);
end

end

%% 
function [grid1, grid2, grid3, isHappy] = ballfunHappiness( vals, pref )
% Check if the function has been resolved. 
    
    r_vals = sum(sum( vals, 2), 3);
    l_vals = sum(sum( vals, 1), 2); 
    t_vals = sum(sum( vals, 1), 3); 
    
    rTech = chebtech2.make( r_vals );
    lTech = trigtech.make( l_vals );
    tTech = trigtech.make( t_vals );
    
    rData.hscale = 1;
    rData.vscale = max( abs( r_vals(:) ) );
    lData.hscale = pi;
    lData.vscale = max( abs( l_vals(:) ) );
    tData.hscale = pi;
    tData.vscale = max( abs( t_vals(:) ) );
    
    resolved_r = happinessCheck(rTech, [], r_vals, rData, pref);
    resolved_l = happinessCheck(lTech, [], l_vals, lData, pref);
    resolved_t = happinessCheck(tTech, [], t_vals, tData, pref);

    isHappy = resolved_r & resolved_l & resolved_t;
    
    [grid1, grid2, grid3] = size(vals);
    if ( ~resolved_r ) 
        grid1 = round( 1.5*size(vals,1) );
    end
    if ( ~resolved_l )
        grid2 = round( 1.5*size(vals,2) );
    end
    if ( ~resolved_t ) 
        grid3 = round( 1.5*size(vals,3) );
    end
    
end

%%
function [isVectorize, op] = vectorCheck(op)
% Check for cases like op = @(x,y,z) x*y^2*z

isVectorize = false;
[xx, yy, zz] = ndgrid([-1,1], [-pi,pi], [-pi,pi]);
try
    A = feval(op, xx, yy, zz);
catch
    throwVectorWarning();
    isVectorize = true;
    return
end

A = feval(op, xx, yy, zz);
if ( any(isinf(A(:) ) ) )
    error('CHEBFUN:BALLFUN:constructor:inf', ...
        'Function returned INF when evaluated');
elseif ( any(isnan(A(:)) ) )
    error('CHEBFUN:BALLFUN:constructor:nan', ...
        'Function returned NaN when evaluated');
end
if ( isscalar(A) )
    op = @(x,y,z) op(x,y,z) + 0*x + 0*y + 0*z;
end
end

%%
function throwVectorWarning()
warning('CHEBFUN:BALLFUN:constructor:vectorize',...
    ['Function did not correctly evaluate on an array.\n', ...
    'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
    'Use the ''vectorize'' flag in the BALLFUN constructor\n', ...
    'call to avoid this warning message.']);
end

%%
function op = str2op(op)
% OP = STR2OP(OP), finds independent variables in a string and returns an
% op handle than can be evaluated.

vars = symvar(op);        % Independent variables
numVars = numel(vars);
if ( numVars == 0 )
    op = @(x,y,z) eval (op);
    
elseif ( numVars == 1 )
    op = eval(['@(' vars{1} ', myVarBeta, myVarGamma)' op]);
    
elseif ( numVars == 2 )
    op = eval(['@(' vars{1} ',' vars{2} ', myVarGamma)' op]);
    
elseif ( numVars == 3 )
    op = eval(['@(' vars{1} ',' vars{2} ',' vars{3} ')' op]);
    
else
    error('CHEBFUN:BALLFUN:constructor:str2op:depvars', ...
        'Too many independent variables in string input.');
end

end