function [Lstar, op] = adjointFormal(L, pref)
%ADJOINTFORMAL   Compute the formal adjoint of a linear differential LINOP.
%   [LSTAR, OP] = ADJOINTFORMAL(L, PREF), where L is a LINOP, returns the formal 
%   adjoint LINOP of L, LSTAR, under the condition that L is a linear 
%   differential operator, i.e. 
%
%   LSTAR*u = a_k*diff(u,k) + ... + a_0*u.
%
%   OP is a function handle that can be used to construct a CHEBOP.
%
% See also ?.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Initialize Lstar and op to have correct dimensions
[nrows ncols] = size(L);
op = mat2cell(zeros(ncols,nrows),ones(1,ncols),ones(1,nrows));
Lstar = linop(op);

%%
% if L.blocks{ii,jj} is a chebfun convert to operator block
for ii = 1:nrows
    for jj = 1:ncols
        if ( strcmp( class(L.blocks{ii,jj}), 'chebfun' ) )
            b = L.blocks{ii,jj};
            L.blocks{ii,jj} = operatorBlock.mult(b,b.domain);
        end
    end
end

%% 
% loop through blocks
op = [];
argstr = '@(x';
for ii = 1:ncols
    for jj = 1:nrows
        % get local block
        B = L.blocks{jj,ii};

        % Get the domain and the value of the highest derivative:
        dom = B.domain;
        dor = B.diffOrder;

        % Get the coefficients:
        coeffs = toCoeff(B, pref);

        % Compute the coefficients of the adjoint:
        adjCoeffs = 0*coeffs;
        for k = 0:dor
            for l = 0:k
                adjCoeffs(dor+1-l) = adjCoeffs{dor+1-l} + ...
                    (-1)^k*nchoosek(k,l)*conj(diff(coeffs{dor+1-k}, k-l));
            end
        end

        % Update the blocks of Lstar using adjoint coefficients
        Lstar.blocks{ii,jj} = 0;
        M = @(f) operatorBlock.mult(f, dom);
        D = @(k) operatorBlock.diff(dom, k);
        for k = 0:dor
            varname = genvarname(['a', int2str(ii), int2str(jj), int2str(k)]);
            eval([varname, '= adjCoeffs{dor+1-k};']);
            Lstar.blocks{ii,jj} = Lstar.blocks{ii,jj} + M(adjCoeffs{dor+1-k}) * D(k);
        end

        % update argstr
        ustr = ['u',int2str(jj)];
        if ( ii == 1 )
            argstr = [argstr,',',ustr];
        end

        % set Bop
        astr = ['a',int2str(ii),int2str(jj)];
        Bop = [astr,'0*',ustr];
        for k = 1:dor
            Bop = [astr,int2str(k),'*diff(',ustr,',',int2str(k),') + ',Bop];
        end
 
        % update op
        if ( jj == 1 )
            op = [op,' ',Bop];
        else
            op = [op,' + ',Bop];
        end

    end

    % update op
    if ( ii < ncols )
        op = [op,';'];
    end

end 

% set op
if ( ncols > 1 )
    op = [argstr,') [',op,'];'];
else
    op = [argstr,') ',op];
end
eval(op);

end
