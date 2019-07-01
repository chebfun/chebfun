function f = constructor( f, op, varargin )
%CONSTRUCTOR   BALLFUN constructor.
%   Given a function OP of three variables, defined on the unit ball, this
%   code represents it as a BALLFUN object. A BALLFUN object is a
%   tensor-product representation of a function that has been "doubled-up"
%   in the radial (r) and polar (theta) variables. This doubled-up function
%   is represented using a Chebyshev-Fourier-Fourier basis.
%
% See also SPHEREFUN, DISKFUN, CHEBFUN2, CHEBFUN3.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 0 )          % BALLFUN( )
    return
end

[op, pref, isVectorized, fixedSize] = parseInputs(op, varargin{:});

% Set preferences:
tech            = pref.tech();
tpref           = tech.techPref;
grid1           = tpref.minSamples;
grid2           = tpref.minSamples;
grid3           = tpref.minSamples;
% Ensure m is odd, n is even, p is even greater than 4
grid1 = grid1 + 1 - mod(grid1,2);
grid2 = grid2 + mod(grid2,2);
grid3 = max(4, grid3 + mod(grid3,2));

maxSample = tpref.maxLength; % maxSample = max grid dimensions.

if ( isa(op, 'ballfun') )     % BALLFUN( BALLFUN )
    f = simplify(op);
    f = fixTheSize(f, fixedSize);
    return
    
elseif ( isa(op, 'double') )   % BALLFUN( DOUBLE )
    f = constructFromDouble(f, op);
    f = fixTheSize(f, fixedSize);
    return
end

isHappy = 0;     % We are currently unresolved.
failure = 0;     % Reached max discretization size without being happy.

while ( ~isHappy && ~failure )
    %% Main loop of the constructor
    vals = evaluate(op, [grid1, grid2, grid3], isVectorized);
    
    % Does the function blow up or evaluate to NaN?:
    vscale = max(abs(vals(:)));
    if ( isinf(vscale) )
        error('CHEBFUN:BALLFUN:constructor:inf', ...
            'Function returned INF when evaluated');
    elseif ( any(isnan(vals(:))) )
        error('CHEBFUN:BALLFUN:constructor:nan', ...
            'Function returned NaN when evaluated');
    end
    
    % If the size of the function is above maxSample then stop.
    if ( max([grid1*grid2, grid2*grid3, grid1*grid3]) > maxSample )
        warning('CHEBFUN:BALLFUN:constructor:dimensions', ...
            'Not well-approximated by a Chebyshev-Fourier-Fourier expansion.');
        failure = 1;
    end
    
    [grid1, grid2, grid3, cutoffs, resolved] = ballfunHappiness( vals );
    isHappy = all(resolved);
    
end

% Evaluate at correct points for the BMC
% Ensure m is odd, n is even, p is even greater than 4
grid1 = cutoffs(1) + 1 - mod(cutoffs(1),2);
grid2 = cutoffs(2) + mod(cutoffs(2),2);
grid3 = max(4, cutoffs(3) + mod(cutoffs(3),2));
[vals, isReal] = evaluate(op, [grid1,grid2,grid3], isVectorized);

% Chop down to correct size: build from coeffs and then cut off

% We are now happy so make a BALLFUN from its values: 
cfs = ballfun.vals2coeffs(vals);

% Simplify: 
if ( resolved(1) )
    cfs = cfs(1:cutoffs(1), :, :); 
end
n = size(cfs,2); 
mid = floor(n/2)+1;
if ( resolved(2) )
    cfs = cfs(:, mid-floor(cutoffs(2)/2):mid+cutoffs(2)-floor(cutoffs(2)/2)-1, :); 
end
p = size(cfs,3); 
mid = floor(p/2)+1;
if ( resolved(3) )
    cfs = cfs(:, :, mid-floor(cutoffs(3)/2):mid+cutoffs(3)-floor(cutoffs(3)/2)-1); 
end

f.coeffs = cfs;
f.isReal = isReal;
f.domain = [0, 1, -pi, pi, 0, pi];
f = fixTheSize(f, fixedSize);
end

%%
function f = constructFromDouble( f, op )
%CONSTRUCTFROMDOUBLE  Construct BALLFUN from matrix of values F = (f_ijk),
% the numbers fijk are used as function values at tensor equally-spaced points 
% in the intrinsic spherical coordinate system, i.e., [0,1]x[-pi,pi]x[0,pi].

if ( numel(op) == 1 )
    f = ballfun(@(x,y,z) op + 0*x);
    return
end

% Double and impose BMC structure
[vals, isReal] = ImposeBMC(op);

% Transform to coefficients
f.coeffs = ballfun.vals2coeffs(vals);
f.isReal = isReal;
f.domain = [0, 1, -pi, pi, 0, pi]; 
f = simplify(f);
end


function [vals, isReal] = evaluate(op, S, isVectorized)
%EVALUATE   Evaluate at a Cheb-Fourier-Fourier grid of size S.
%  EVALUATE(g, S) returns the S(1)xS(2)xS(3) values of g at a
%  Chebyshev-Fourier-Fourier grid for the function
%  g(r, lambda, theta).

% Convert a handle_function to a ballfun function
m = S(1);
n = S(2);
p = S(3);

% Evaluation points (assuming m odd, n even and p even > 4)
r = chebpts(m);
lam = [pi*trigpts(n); pi];
th = [pi*trigpts(p); pi];

% Add dependence on r, lambda and theta
f1 = @(r,lam,th)op(r,lam,th) + 0*r + 0*lam + 0*th;

% Build the grid of evaluation points
[rrg, llg, ttg] = ndgrid(r(floor(m/2)+1:m), lam(1:n/2+1), th(p/2+1:p+1));
[rrh, llh, tth] = ndgrid(r(floor(m/2)+1:m), lam(n/2+1:end), th(p/2+1:p+1));

if ( ~isVectorized )
    % g : evaluation at [0,1] x [-pi,0] x [0,pi]
    g = feval(f1, rrg, llg, ttg);
    % h : evaluation at [0,1] x [0,pi] x [0,pi]
    h = feval(f1, rrh, llh, tth);
else
    % If vectorize flag is turned on, then FOR loop:
    % g and h have the same size
    g = zeros(size(rrg));
    h = zeros(size(rrh)); 
    for i1 = 1:size(g,1)
        for j1 = 1:size(g,2)
            for k1 = 1:size(g,3)
                % g : evaluation at [0,1] x [-pi,0] x [0,pi]
                g(i1,j1,k1) = f1( rrg(i1,j1,k1), llg(i1,j1,k1), ttg(i1,j1,k1) );
                % h : evaluation at [0,1] x [0,pi] x [0,pi]
                h(i1,j1,k1) = f1( rrh(i1,j1,k1), llh(i1,j1,k1), tth(i1,j1,k1) ); 
            end
        end
    end
end

% Double the function in r and theta
[vals, isReal] = ImposeBMC(g,h);
end

function [vals, isReal] = ImposeBMC(f,h)
% Take a tensor of values on [0,1]x[-pi,pi[,[0,pi] and double it in the
% direction r and theta, then impose the BMCIII structure

if nargin == 1
    % Get the size
    [m,n,p] = size(f);
    
    if ( mod(n,2) ~= 0 )
        error('CHEBFUN:BALLFUN:constructor:sizevalues', ... 
        'When constructing from values the number of columns must be even.');
    end
    
    if ( p < 2 )
        error('CHEBFUN:BALLFUN:constructor:sizevalues', ... 
        'When constructing from values the number of tubes must be greater than 2.');
    end
  
    % Define the functions g on [0,1]x[-pi,0]x[0,pi] and h on [0,1]x[0,pi]x[0,pi]
    g = f(:,1:n/2+1,:);
    h = f(:,[n/2+1:end 1],:); 
else
    % Get the size
    [m,n,p] = size(f);
    
    g = f;
    % Doubled size in theta
    n = 2*n-2;
end

% Doubled size in r
m = 2*m-1;

% Doubled size in theta
p = 2*p-2;

vals = zeros(m,n,p);

%% Impose BMC-III Structure
% f(0,:,:) = constant
% Count lambda = pi only once
% Compute the mean of f evaluated at r = 0
g0 = g(1,:,:); h0 = h(1,2:end,:);
m_zeroR = mean([mean(g0(:)), mean(h0(:))]);
g(1,:,:) = m_zeroR;
h(1,:,:) = m_zeroR;

% f(r,:,0) = constant
m_zeroT = mean([mean(g(:,:,1),2),mean(h(:,2:end,1),2)],2);
g(:,:,1) = repmat(m_zeroT,1,size(g,2));
h(:,:,1) = repmat(m_zeroT,1,size(g,2));

% f(r,:,pi) = constant
m_piT = mean([mean(g(:,:,end),2),mean(h(:,2:end,end),2)],2);
g(:,:,end) = repmat(m_piT,1,size(g,2));
h(:,:,end) = repmat(m_piT,1,size(g,2));

%% Flip g and h on the radial direction
flip1g = flip(g(1+mod(m,2):end,:,:), 1);
flip1h = flip(h(1+mod(m,2):end,:,:), 1);

%% Fill in vals
% [0,1] x [-pi,0] x [0,pi[
vals(floor(m/2)+1:m, 1:n/2+1, floor((p+1)/2)+1:p) = g(:,:,1:end-1);
% [0,1] x [0,pi[ x [0,pi[
vals(floor(m/2)+1:m, n/2+1:n, floor((p+1)/2)+1:p) = h(:,1:end-1,1:end-1);
% [-1,0[ x [-pi,0] x [0,pi[
vals(1:floor(m/2), 1:n/2+1, floor((p+1)/2)+1:p) = flip(flip1h(:,:,2:end),3);
% [-1,0[ x [0,pi[ x [0,pi[
vals(1:floor(m/2), n/2+1:n, floor((p+1)/2)+1:p) = flip(flip1g(:,1:end-1,2:end),3);
% [0,1] x [-pi,0] x [-pi,0]
vals(floor(m/2)+1:m, 1:n/2+1, 1:floor((p+1)/2)) = flip(h(:,:,2:end),3);
% [0,1] x [0,pi[ x [-pi,0]
vals(floor(m/2)+1:m, n/2+1:n, 1:floor((p+1)/2)) = flip(g(:,1:end-1,2:end),3);
% [-1,0[ x [0,pi[ x [-pi,0]
vals(1:floor(m/2), n/2+1:n, 1:floor((p+1)/2)) = flip1h(:,1:end-1,1:end-1);
% [-1,0[ x [-pi,0] x [-pi,0]
vals(1:floor(m/2), 1:n/2+1, 1:floor((p+1)/2)) = flip1g(:,:,1:end-1);

% Check if the function is real
if norm(imag(vals(:))) < 10^5*eps
   isReal = 1;
   vals = real(vals);
else
    isReal = 0;
end
end

%%
function f = fixTheSize(f, fixedSize)
% Fix the size of f to the one provided by fixedSize
if ~isempty(fixedSize)
   f.coeffs = coeffs3(f, fixedSize(1), fixedSize(2), fixedSize(3));
end
end

%%
function [op, pref, isVectorized, fixedSize] = parseInputs(op, varargin)
% Parse user inputs to BALLFUN.

isVectorized = 0;
isCoeffs = 0;
fixedSize = [];
pref = chebfunpref();

% Preferences structure given?
isPref = find(cellfun(@(p) isa(p, 'chebfunpref'), varargin));
if ( any(isPref) )
    pref = varargin{isPref};
    varargin(isPref) = [];
end

if ( isa(op, 'char') )     % BALLFUN( CHAR )
    op = str2op(op);
end

if ( isa(op, 'function_handle') )
    % Check for polar coords
    ispolar = any(find(strcmp(varargin,'polar'))) || any(find(strcmp(varargin,'spherical')));
    if ( ~ispolar )
        x = @(r,lam,th)r.*sin(th).*cos(lam);
        y = @(r,lam,th)r.*sin(th).*sin(lam);
        z = @(r,lam,th)r.*cos(th);
        op = @(r,lam,th) op(x(r,lam,th), y(r,lam,th), z(r,lam,th));
    end
end

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    elseif any(strcmpi(varargin{k}, {'vectorize', 'vectorise'}))
        isVectorized = true;
    elseif strcmpi(varargin{k}, 'coeffs')
        isCoeffs = 1;
    elseif isa(varargin{k}, 'double') && length(varargin{k}) == 3
        fixedSize = varargin{k};
    end
end

if ( isCoeffs )
    % We have all we need:
    cfs = op;
    % Test if the matrix is sparse
    if issparse(cfs)
        cfs = full(cfs);
    end
    op = ballfun();
    op.coeffs = cfs;
    
    % Check if the function is real
    [~,n,p] = size(cfs);
    % Check if f = conj(f) on the coeffs
    CheckReal = cfs(:,2-mod(n,2):floor(n/2)+1,2-mod(p,2):floor(p/2)+1)-conj(cfs(:,end:-1:floor(n/2)+1,end:-1:floor(p/2)+1));
    bReal = norm(CheckReal(:),inf) < 10^7*eps;
    % Additional check if n or p is even
    if mod(n,2) == 0
        CheckReal = cfs(:,1,:);
        bReal = bReal && ( norm(CheckReal(:),inf) < 10^7*eps ); 
    end
    if mod(p,2) == 0
        CheckReal = cfs(:,:,1);
        bReal = bReal && ( norm(CheckReal(:),inf) < 10^7*eps );
    end
    op.isReal = bReal;
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~isVectorized && ~isnumeric(op) && ~isCoeffs) % another check
    [isVectorized, op] = vectorCheck(op);
end
end

%% 
function [grid1, grid2, grid3, cutoffs, resolved] = ballfunHappiness( vals )
% Check if the function has been resolved. 
    
    vscl = max(1, max( abs( vals(:) ) )); 
    cfs = ballfun.vals2coeffs( vals ); 

    r_cfs = max(max( abs(cfs), [], 2), [], 3);
    l_cfs = max(max( abs(cfs), [], 1), [], 3);
    l_cfs = l_cfs(:);
    t_cfs = max(max( abs(cfs), [], 1), [], 2);
    t_cfs = t_cfs(:);
    
    rTech = chebtech2.make( {'',r_cfs} );
    lTech = trigtech.make( {'',l_cfs} );
    tTech = trigtech.make( {'',t_cfs} );
    
    rvals = rTech.coeffs2vals(rTech.coeffs);
    rdata.vscale = vscl; 
    rdata.hscale = 1;
    lvals = lTech.coeffs2vals(lTech.coeffs);
    ldata.vscale = vscl; 
    ldata.hscale = 1;
    tvals = tTech.coeffs2vals(tTech.coeffs);
    tdata.vscale = vscl; 
    tdata.hscale = 1;
    
    % Check happiness along each slice: 
    [resolved_r, cutoff_r] = happinessCheck(rTech, [], rvals, rdata);
    [resolved_l, cutoff_l] = happinessCheck(lTech, [], lvals, ldata);
    [resolved_t, cutoff_t] = happinessCheck(tTech, [], tvals, tdata);
    
    [grid1, grid2, grid3] = size(vals);
    if ( ~resolved_r ) 
        grid1 = round( 1.5*size(vals,1) );
        grid1 = grid1 + 1 - mod(grid1,2);
    end
    if ( ~resolved_l )
        grid2 = round( 1.5*size(vals,2) );
        grid2 = grid2 + mod(grid2,2);
    end
    if ( ~resolved_t ) 
        grid3 = round( 1.5*size(vals,3) );
        grid3 = max(4, grid3 + mod(grid3,2));
    end
    cutoffs = [cutoff_r, cutoff_l, cutoff_t];
    resolved = [resolved_r, resolved_l, resolved_t];
end

%%
function [isVectorized, op] = vectorCheck(op)
% Check for cases like op = @(x,y,z) x*y^2*z

isVectorized = false;
[xx, yy, zz] = ndgrid([-1,1], [-pi,pi], [-pi,pi]);
try
    A = feval(op, xx, yy, zz);
catch
    throwVectorWarning();
    isVectorized = true;
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