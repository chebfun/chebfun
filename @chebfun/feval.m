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

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Add support for multi-rowed impulses and deltas and update the
% preceding text accordingly.

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
    if ( any(strcmpi(x, {'left', 'start' ,'-'})) )
        out = get(F, 'lval');
    elseif ( any(strcmpi(x, {'right', 'end', '+'})) )
        out = get(F, 'rval');
    else
        error('CHEBFUN:feval:strInput', 'Unknown input argument "%s".', x);
    end
    return
end

%% DEAL WITH QUASIMATRICES:
out = cell(1, numel(F));
if ( F(1).isTransposed )
    out = out.';
end
for k = 1:numel(F)
    out{k} = columnFeval(F(k), x, varargin{:});
end
out = cell2mat(out);
        
end

function out = columnFeval(f, x, varargin)

%% INITIALISE:
% Un-transpose f if necessary so that the call to size() below returns the
% right thing.  (We'll compensate for this further down.)
wasTransposed = f.isTransposed;
if ( f.isTransposed )
    f = transpose(f);
end

% Reshape x to be a column vector.
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

% Initialise output:
numCols = size(f, 2);
numFuns = numel(f.funs);
out = zeros(size(x, 1), numCols);
funs = f.funs;
dom = f.domain;

%% LEFT AND RIGHT LIMITS:
% Deal with feval(f, x, 'left') and feval(f, x, 'right'):
if ( nargin > 2 )
    lr = varargin{1};
    lrFlag = strcmpi(lr, {'left', 'right', '-', '+'});
    if ( ~any(lrFlag) )
        if ( ischar(lr) )
            error('CHEBFUN:feval:leftrightchar',...
                'Unknown input argument "%s".', lr);
        else
            error('CHEBFUN:feval:leftright', 'Unknown input argument.');
        end
    end
    % We deal with this by reassigning imps to be left/right values.
    if ( lrFlag(1) || lrFlag(3) ) % left
        f.impulses(:,:,2:end) = []; % Level 2 imps are not needed here
        for j = 1:numFuns
            f.impulses(j+1,:,1) = get(funs{j}, 'rval');
        end
    elseif ( lrFlag(2) || lrFlag(4) ) % right
        f.impulses(:,:,2:end) = []; % Level 2 imps are not needed here
        for j = 1:numFuns
            f.impulses(j,:,1) = get(funs{j}, 'lval');
        end
    end
end

%% VALUES FROM FUNS:

% Points to the left of the domain:
xReal = real(x);
I = xReal < dom(1);
if ( any(I(:)) )
    out(I,:) = feval(funs{1}, x(I));
end

% Points within the domain:
for k = 1:numFuns
    I = ( xReal >= dom(k) ) & ( xReal < dom(k+1) );
    if ( any(I(:)) )
        % Evaluate the appropriate fun:
        out(I,:) = feval(funs{k}, x(I));
    end
end

% Points to the right of the domain:
I = ( xReal >= dom(end) );
if ( any(I(:)) )
    out(I,:) =  feval(funs{end}, x(I));
end

%% IMPULSES:
% If the evaluation point corresponds to a breakpoint, we get the value from
% imps. If there is only one row, the value is given by the corresponding entry
% in that row. If the second row is nonzero the value is -inf or inf
% corresponding to the sign of the entry in the 2nd row. If the entry in the
% corresponding 3rd or higher rows is nonzero, we return NaN.

higherImpulses = f.impulses(:,:,2:end);
% Only one row:
if ( ~any(higherImpulses(:)) )
    % Loop over the FUNs:
    for k = 1:numFuns + 1
        index = x == dom(k);
        if ( any(index) )
            imps = repmat(f.impulses(k,:,1), sum(index), 1);
            out(index,:) = imps;
        end
    end
    
else
    % [TODO]: Higher-order impulses.
end

%% RESHAPE FOR OUTPUT:
% Reshape fx, which is a column vector or horizontal concatenation of column
% vectors, to be of the appropriate size, and handle transposition.
sizefx = sizex;
sizefx(2) = numCols*sizex(2);
if ( ndimsx == 2 )
    % If x was just a matrix or vector, the reshape is straightforward.
    out = reshape(out, sizefx);

    if ( wasTransposed )
        out = out.';
    end
else
    % If x had more than two dimensions, we have to be more careful.  The
    % cell2mat(mat2cell(...).') effects a transpose which keeps certain
    % columnar blocks of the fx matrix intact, putting the entries in the
    % correct order for reshape().
    blockLength = sizex(1)*sizex(2);
    blocksPerCol = prod(sizex(3:end));
    out = reshape(cell2mat(mat2cell(out, blockLength*ones(1, blocksPerCol), ...
        ones(1, numCols)).'), sizefx);

    % We define "transposition" in this case to mean the switching of the first
    % two dimensions.  [TODO]:  Is this the right thing to do?
    if ( wasTransposed )
        out = permute(out, [2 1 3:ndimsx]);
    end
end

end
