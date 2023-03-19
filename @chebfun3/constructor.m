function f = constructor(f, op, varargin)
%CONSTRUCTOR   CHEBFUN3 constructor.
%   Given a function OP of three variables, this code represents it as a
%   CHEBFUN3 object. A CHEBFUN3 object is a low-rank representation
%   expressing a function as a trilinear product of a discrete core tensor
%   and three quasimatrices consisting of univariate functions.
%
%
%   Since March 2023, chebfun3f is used by default to construct CHEBFUN3
%   objects. The legacy constructor chebfun3classic can be accessed by
%   providing the 'classic' flag as additional input argument.
%
% See also CHEBFUN2, CHEBFUN3T and CHEBFUN3V.

% Copyright 2023 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse the inputs:
[op, dom, pref, fiberDim, vectorize, isEqui, fixedRank] = ...
    parseInputs(op, varargin{:});

% The 'equi' flag can be used only with numeric data:
if ( isEqui && ~isa(op, 'double') )
    error('CHEBFUN:CHEBFUN3:constructor:equi', ...
        'The EQUI flag is valid only when constructing from numeric data');
end

if ( isa(op, 'chebfun3') )     % CHEBFUN3( CHEBFUN3 )
    f = op;
    return
elseif ( isa(op, 'double') )   % CHEBFUN3( DOUBLE )
    f = chebfun3.chebfun3double(f, op, dom, pref, isEqui);
elseif ( strcmpi(pref.cheb3Prefs.constructor, 'classic') )
    f = chebfun3.chebfun3classic(f, op, pref, dom, vectorize, fiberDim);
else
    f = chebfun3.chebfun3f(f, op, pref, dom, vectorize);
end

if ( fixedRank )
    % Simplify the rank if requested.
    f = fixTheRank(f , fixedRank);
end
end
%% End of constructor

%%
function [op, dom, pref, fiberDim, vectorize, isEqui, ...
    fixedRank] = parseInputs(op, varargin)
vectorize = 0;
isEqui = 0;
isCoeffs = 0;
fixedRank = 0;
pref = chebfunpref();
fiberDim = [];
dom = [-1 1 -1 1 -1 1];

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
    if any(strcmpi(varargin{k}, {'trig', 'periodic'}))
        pref.tech = @trigtech;
    elseif strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    elseif strcmpi(varargin{k}, 'classic')
        pref.cheb3Prefs.constructor = 'classic';
    elseif strcmpi(varargin{k}, 'chebfun3f')
        pref.cheb3Prefs.constructor = 'chebfun3f';
    elseif strcmpi(varargin{k}, 'rank') % rank is specified.
        fixedRank = varargin{k+1};
    elseif ( isnumeric(varargin{k}) )
        if ( numel(varargin{k}) == 6 ) % domain is specified.
            dom = varargin{k};
        elseif ( numel(varargin{k}) == 3 ) % length is specified.
            if ( k > 1 && strcmpi(varargin{k-1}, 'rank') )
                % Rank is specified, not length. Don't confuse them.
                continue
            else
                % Interpret this as the user wants a fixed degree chebfun3
                % on the domain DOM.
                len = varargin{k};
                [xx, yy, zz] = chebpts3(len(1), len(2), len(3), dom);
                op = op(xx, yy, zz);
            end
        end
    elseif strcmpi(varargin{k}, {'fiberDim'})
        fiberDim = varargin{k+1};
    elseif any(strcmpi(varargin{k}, {'vectorize', 'vectorise'}))
        vectorize = true;
    elseif strcmpi(varargin{k}, 'coeffs')
        isCoeffs = 1;
    elseif strcmpi(varargin{k}, 'equi')
        isEqui = 1;
    end
end

if ( isCoeffs )
    op = chebfun3.coeffs2vals(op);
end

% If the vectorize flag is off, do we need to give user a warning?
if ( ~vectorize && ~isnumeric(op) ) % another check
    [vectorize, op] = vectorCheck(op, dom);
end

end

%%
function [vectorize, op] = vectorCheck(op, dom)
% Check for cases like op = @(x,y,z) x*y^2*z

vectorize = false;
[xx, yy, zz] = ndgrid(dom(1:2), dom(3:4), dom(5:6));
try
    A = feval(op, xx, yy, zz);
catch
    throwVectorWarning();
    vectorize = true;
    return
end

A = feval(op, xx, yy, zz);
if ( any(isinf(A(:) ) ) )
    error('CHEBFUN:CHEBFUN3:constructor:inf', ...
        'Function returned INF when evaluated');
elseif ( any(isnan(A(:)) ) )
    error('CHEBFUN:CHEBFUN3:constructor:nan', ...
        'Function returned NaN when evaluated');
end

if ( isscalar(A) )
    op = @(x,y,z) op(x,y,z) + 0*x + 0*y + 0*z;
end

end

%%
function throwVectorWarning()
warning('CHEBFUN:CHEBFUN3:constructor:vectorize',...
    ['Function did not correctly evaluate on an array.\n', ...
    'Turning on the ''vectorize'' flag. Did you intend this?\n', ...
    'Use the ''vectorize'' flag in the CHEBFUN3 constructor\n', ...
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
    error('CHEBFUN:CHEBFUN3:constructor:str2op:depvars', ...
        'Too many independent variables in string input.');
end

end


%%
function f = fixTheRank(f , fixedRank)
% Fix the rank of a CHEBFUN3. Used for calls to constructor with specified
% rank.

if ( any(fixedRank) < 0 )
    error('CHEBFUN:CHEBFUN3:constructor:fixTheRank:negative', ...
        'Ranks should all be nonnegative.')
elseif ( all(fixedRank) )
    [r1, r2, r3] = rank(f);
    t1 = fixedRank(1);
    t2 = fixedRank(2);
    t3 = fixedRank(3);
    
    % What to do with cols?
    if ( r1 > t1 )
        % Truncate cols:
        f.cols = f.cols(:, 1:t1);
        f.core = f.core(1:t1, :, :);
        r1 = t1; % New size of core
    elseif ( r1 < t1 )
        % Pad cols with approprate number of zero cols:
        zCols = chebfun(0, f.cols.domain);
        for jj = r1 : t1 - 1
            f.cols = [f.cols zCols];
        end
        % Pad mode 1 of the core tensor with zeros:
        tempCore = zeros(t1, r2, r3);
        tempCore(1:r1, :, :) = f.core;
        f.core = tempCore;
        r1 = t1; % New size of core
    end
    
    % What to do with rows?
    if ( r2 > t2 )
        % Truncate rows:
        f.rows = f.rows(:,1:t2);
        f.core = f.core(:, 1:t2, :);
        r2 = t2; % New size of core
    elseif ( r2 < t2 )
        % Pad rows with approprate number of zero rows:
        zRows = chebfun(0, f.rows.domain);
        for jj = r2 : t2 - 1
            f.rows = [f.rows zRows];
        end
        % Pad mode 2 of the core tensor with zeros:
        tempCore = zeros(r1, t2, r3);
        tempCore(:, 1:r2, :) = f.core;
        f.core = tempCore;
        r2 = t2; % New size of core
    end
    
    % What to do with tubes?
    if ( r3 > t3 )
        % Truncate tubes:
        f.tubes = f.tubes(:, 1:t3);
        f.core = f.core(:, :, 1:t3);
    elseif ( r3 < t3 )
        % Pad tubes with approprate number of zero tubes:
        zTubes = chebfun(0, f.tubes.domain);
        for jj = r3 : t3 - 1
            f.tubes = [f.tubes zTubes];
        end
        % Pad mode 3 of the core tensor with zeros:
        tempCore = zeros(r1, r2, t3);
        tempCore(:, :, 1:r3) = f.core;
        f.core = tempCore;
    end
end

end



