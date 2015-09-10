function out = feval(F, x, varargin)
%FEVAL   Evaluate a CHEBFUN.
%   FEVAL(F, X) evaluates a CHEBFUN F at the points in X.  If F is a quasimatrix
%   with columns F1, ..., FN, then the result will be [F1(X), ..., FN(X)], the
%   horizontal concatenation of the results of evaluating each column at the
%   points in X.
%
%   FEVAL(F, 'left'), FEVAL(F, 'start'), and FEVAL(F, '-') return the value of F
%   at the left endpoint of its domain.  FEVAL(F, 'right'), FEVAL(F, 'end'), and
%   FEVAL(F, '+') do the same for the right endpoint.
%
%   FEVAL(F, X, 'left') and FEVAL(F, X, '-') evaluate F at the points in X,
%   using left-hand limits to evaluate F at any breakpoints. FEVAL(F, X,
%   'right') and FEVAL(F, X, '+') do the same but using right-hand limits.
%
%   F(X), F('left'), F(X, 'left'), etc, are equivalent syntaxes. 
%
%   Example:
%     f = chebfun(@(x) 1./(1 + 25*x.^2));
%     y = feval(f, linspace(-1, 1, 100).');
%
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If F or x are empty, there's nothing to do.
if ( isempty(F) )
    out = [];
    return
elseif ( isempty(x) )
    % Return empty matrix with dimensions of the appropriate size.
    out = zeros(size(x));
    return
end

if ( isa(F, 'function_handle') )
    out = F(x, varargin{:});
    return
end

%% LEFT / RIGHT VALUES:
% Support for feval(f, 'left') and feval(f, 'end'), etc.
if ( ischar(x) )
    dom = F.domain;
    if ( any(strcmpi(x, {'left', 'start' ,'-'})) )
        out = feval(F, dom(1));
    elseif ( any(strcmpi(x, {'right', 'end', '+'})) )
        out = feval(F, dom(end));
    else
        error('CHEBFUN:CHEBFUN:feval:strInput', ...
            'Unknown input argument "%s".', x);
    end
    return
end

%% DEAL WITH QUASIMATRICES:
out = cell(1, numel(F));
for k = 1:numel(F)
    out{k} = columnFeval(F(k), x, varargin{:});
end
out = cell2mat(out);
if ( F(1).isTransposed )
    % We got a passed a row CHEBFUN. If X had more than two dimensions, we can't
    % simply transpose the output from above, instead, we need to use permute.
    ndimsx = ndims(x);
    if ( ndimsx <= 2 )
        out = out.';
    else
        % We define "transposition" in this case to mean the switching of the
        % first two dimensions.
        out = permute(out, [2 1 3:ndimsx]);
    end
end
        
end

function out = columnFeval(f, x, varargin)

%% INITIALISE:

% Reshape x to be a column vector.
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

% Initialise output:
numCols = numColumns(f);
numFuns = numel(f.funs);

funs = f.funs;
dom = f.domain;

%% LEFT AND RIGHT LIMITS:
% Deal with feval(f, x, 'left') and feval(f, x, 'right'):
lrFlag = zeros(1, 4);
if ( nargin > 2 )
    lr = varargin{1};
    lrFlag = strcmpi(lr, {'left', 'right', '-', '+'});
    if ( ~any(lrFlag) )
        if ( ischar(lr) )
            error('CHEBFUN:CHEBFUN:feval:leftRightChar',...
                'Unknown input argument "%s".', lr);
        else
            error('CHEBFUN:CHEBFUN:feval:leftRight', 'Unknown input argument.');
        end
    end
end

%% VALUES FROM FUNS:

if ( numFuns == 1 )
    
    % Things are simple when there is only a single FUN:
    out = feval(funs{1}, x(:), varargin{:});
    
else
    
    % For multiple FUNs we must determine which FUN corresponds to each x.
    
    % Initialise output matrix:
    out = zeros(numel(x), numCols);
    
    % Replace the first and last domain entries with +/-inf. (Since we want to
    % use FUN{1} if real(x) < dom(1) and FUN{end} if real(x) > dom(end)).
    domInf = [-inf, dom(2:end-1), inf];
    
    % Loop over each fun. If  real(x) is in [dom(k) dom(k+1)] then use FUN{k}.
    xReal = real(x);
    for k = 1:numFuns
        I = ( xReal >= domInf(k) ) & ( xReal < domInf(k+1) );
        if ( any(I(:)) )
            % Evaluate the appropriate fun:
            out(I,:) = feval(funs{k}, x(I), varargin{:});
        end
    end
    
end

%% POINTVALUES:
% If the evaluation point corresponds to a breakpoint, we get the value from
% pointValues. 

% Loop over the FUNs:
for k = 1:(numFuns + 1)
    index = ( x == dom(k) );
    if ( any(index) )
        % If a left or right flag has been passed, we reassign pointValues 
        % to be left/right values.
        if ( lrFlag(1) || lrFlag(3) ) % left
            for j = 1:numFuns
                f.pointValues(j+1,:) = get(funs{j}, 'rval');
            end
        elseif ( lrFlag(2) || lrFlag(4) ) % right
            for j = 1:numFuns
                f.pointValues(j,:) = get(funs{j}, 'lval');
            end
        end
        pointValues = repmat(f.pointValues(k,:), sum(index), 1);
        out(index,:) = pointValues;
    end
end

%% RESHAPE FOR OUTPUT:
% Reshape fx, which is a column vector or horizontal concatenation of column
% vectors, to be of the appropriate size, and handle transposition.
sizefx = sizex;
sizefx(2) = numCols*sizex(2);
if ( ndimsx == 2 )
    % If x was just a matrix or vector, the reshape is straightforward.
    out = reshape(out, sizefx);
else
    % If x had more than two dimensions, we have to be more careful.  The
    % cell2mat(mat2cell(...).') effects a transpose which keeps certain
    % columnar blocks of the fx matrix intact, putting the entries in the
    % correct order for reshape().
    blockLength = sizex(1)*sizex(2);
    blocksPerCol = prod(sizex(3:end));
    out = reshape(cell2mat(mat2cell(out, blockLength*ones(1, blocksPerCol), ...
        ones(1, numCols)).'), sizefx);
end

end
