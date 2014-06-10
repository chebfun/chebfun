function [M, S] = instantiate(disc)
%INSTANTIATE   Convert a ULTRAS discretization to discrete form.
%   M = INSTANTIATE(DISC) converts each item DISC.SOURCE to discrete form
%   using the information in discretization DISC. The result M is return a cell
%   array if DISC.SOURCE has more than one component.
%
%   [M, S] = INSTANTIATE(DISC) retusn a second output, S, which is a cell array
%   containing the dscrete form of the ultraS conversion operator for each block
%   of DISC.SOURCE.
%
%   DISC.SOURCE may be one or a cell array of:
%      linBlock (becomes a matrix)
%      chebfun (becomes a vector)
%      numeric (not changed)

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

data = disc.source;
if ( isa(data, 'chebmatrix') )
    data = data.blocks;
end

if ( iscell(data) )
    M = cell(size(data));
    S = cell(size(data));
    for j = 1:size(data, 1)
        for k = 1:size(data, 2)
            discJK = extractBlock(disc, j, k);
            [M{j,k}, S{j,k}] = instantiate(discJK);
        end
    end
    return
else
    [M, S] = instantiateOne(disc, data);
end

end

function [M, S] = instantiateOne(disc, item)
% Instantiate one block of data.

if ( isa(item, 'operatorBlock') )
    % Convert a square block
    
    if ( ~isempty(disc.coeffs) )
        % Coefficients of the block are available, convert to a diffmat.
        [M, S] = quasi2USdiffmat(disc);
    else
        error('CHEBFUN:ultraS:fail', ...
            'ultraS cannot represent this operator. Suggest you use colloc2.')
    end
    
elseif ( isa(item, 'functionalBlock') )
    % Convert a row block.
    
    % Developer note: In general we can't represent functional
    % blocks via coeffs. To get around this we instantiate a
    % COLLOC2 discretization and convert it to coefficient space
    % using COEFFS2VALS(). (Note it's COEFFS2VALS() rather than
    % VALS2COEFFS() because it's a right-multiply (I think..).)
    

    
    % For convenience:
    dim = disc.dimension;
    dom = disc.domain;   
    
%     [p, loc] = recoverCoeffs( linop(item) );
%     
%     if ( item.diffOrder > 0 && all(~isinf(loc{:})) )
%         % Try simple thing first: 
%         cumsumDim = [0, cumsum(dim)];
%         tmp = cell(1, numel(dom)-1);
%         for l = 1:numel(tmp)
%             bcn = cumsumDim(l) + dim(l);
%             tmp{l} = p{l}*chebValues(0:(numel(p)-1), bcn, loc{l}).';
%         end
%         M = cell2mat(tmp);
%     else
        % Create a colloc2 discretization:
        collocDisc = colloc2(item, dim, dom);
        M = matrix(collocDisc);
    
        % Convert from colloc-space to coeff-space using COEFFS2VALS.
        cumsumDim = [0, cumsum(dim)];
        tmp = cell(1, numel(dom)-1);
        for l = 1:numel(tmp)
            Ml = M(cumsumDim(l) + (1:dim(l)));
            tmp{l} = flipud(chebtech2.coeffs2vals(Ml.')).';
        end
    
        M = cell2mat(tmp);
%     end
    S = zeros(size(M));
    
elseif ( isa(item, 'chebfun') )
    % Block is a CHEBFUN. Convert to value space.
    
    M = toValues(disc, item);
    if ( item.isTransposed )
        M = M.';
    end
    S = zeros(size(M));
    
elseif ( isnumeric(item) )
    % Block is numeric, don't need to do much.
    
    M = item;
    S = 1;
    
else
    
    error('CHEBFUN:ultraS:instantiate:inputType', ...
        'Unrecognized item type.')
    
end

end



function [p, loc] = recoverCoeffs( L )
%RECOVERCOEFFS  Recover coefficient functions of a linear operator
% P = RECOVERCOEFFS(L) returns, for a linear operator L, a chebfun
% quasimatrix P such that
%         Lu = P(:,1)*u + P(:,2)*u' + P(:,3)*u" + ... P(:,M+1)*u^(M),
% where M is the difforder of the operator. If L is not linear, an error is
% thrown.
%
% For a block operator L, i.e., one defining a system of equations
%         Lu = [L_{1,1} L_{1,2} ... L_{1,S}] [ u_1 ]
%              [L_{2,1} L_{1,2} ... L_{1,S}] [ u_2 ]
%              [  ...     ...   ...   ...  ] [ ... ]
%              [L_{R,1} L_{R,2} ... L_{R,S}] [ u_S ],
% P will be the RxS cell array such that P{J,K} = RECOVERCOEFFS(L_{J,K}).
%
% [P L] = RECOVERCOEFFS(L) returns also the linop L, which can be useful if
% the input was a linear chebop.
%
% Example 1:
%  [L x] = chebop(@(x,u) 0.5*diff(u,2) - sin(x).*diff(u) + x.*u);
%  p = recoverCoeffs(L)
%  norm(p - [x -sin(x) 0.5])
%
% Example 2:
%  [L x] = chebop(@(x,u) diff(sin(x).*(diff(cos(x).*u))),[-pi pi]);
%  p = recoverCoeffs(L)
%  norm(p - [-sin(2*x) 1-3*sin(x).^2 sin(2*x)/2])
%
% Example 3:
%  L = chebop(@(x,u,v) [diff(u,2), 0.5*diff(v)+exp(x)]);
%  p = recoverCoeffs(L)
%  norm([p{:}] - [0 0 1 0 0 0 .5])

% Convert to linop if input is a chebop. (But don't overwrite input as it's
% more efficient to evaluate the chebop .op than the linearised .oparray!)
if isa(L,'chebop'), L2 = linop(L); else L2 = L; end

% Initialise
s = size( L2 );                  % Determine the size of the system
m = L2.diffOrder;                %  and the difforder
x = chebfun('x',L2.domain,2);    % Construct linear function on the domain
x0 = chebfun(0,L2.domain);       %  and the zero function
p = cell(s);                     % Initialise output
p0 = L*repmat(x0,1,s(2));        % Compute non-autonomous component

% The main routine
for hh = 1:s(2)                 % Loop over each of the dependant variables
    x0l = repmat(x0,1,hh-1);    % Set dep vars to the left to zero
    x0r = repmat(x0,1,s(2)-hh); % Set dep vars to the right to zero
    p1 = L*[x0l 1+0*x x0r];     % Evaluate all equations for [0 ... 1 ...0]
    p1 = p1 - p0;               % Subtract non-autonomous compnent
    for ll = 1:s(1)             % Loop over equations and assign
        p{ll,hh} = p1{ll};
    end
    xk = x;                            % Update indep var to x
    for kk = 1:max(m(:,hh))            % Loop over each x^k
        tmp = L*[x0l xk(:,kk) x0r]-p0; % Evaluate for u = [0 ... x^k ... 0]
        for ll = 1:s(1)                % Loop over each equation
            if kk > m(ll,hh), continue, end % No coeffs of this order here
            p{ll,hh}(:,kk+1) = tmp{ll};   % Assign the ll-th equation
            for jj = 1:kk              % Extract the lower-order terms
                p{ll,hh}(:,kk+1) = p{ll,hh}(:,kk+1) - p{ll,hh}(:,kk+1-jj).*xk(:,jj);
                p{ll,hh}(:,kk+1) = simplify(p{ll,hh}(:,kk+1)); % Simplify
            end
        end
        xk = [xk x.*xk(:,end)/(kk+1)]; % Update indep var to x^k/k!
    end
    tmp = feval(L2, x.^(numel(p))/factorial(numel(p)));
    loctemp =  tmp{:} / p{end,hh};
    dom = L2.domain; 
    affine = @(s) 2*(s - dom(1))/(dom(2)-dom(1)) - 1;
    loc{hh} = affine(loctemp);
end

end

function val = chebValues(k, n, pos)
% CHEBBALUES. Return the values of ChebT and its derivatives. 
if ( k == 0 )
    val = cos((0:(n-1))'*acos(pos)); %pos.^((0:n-1).');
else
    [ll, kk] = meshgrid((0:n-1),(0:k-1));
    val = (pos).^((1:n).').*prod( (ll.^2 - kk.^2)./(2*kk+1), 1 ).';
end
end

