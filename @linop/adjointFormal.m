function [Lstar, op] = adjointFormal(L, pref)
%ADJOINTFORMAL   Compute the formal adjoint of a LINOP.
%   [LSTAR, OP] = ADJOINTFORMAL(L, PREF), where L is a LINOP, returns the
%   formal adjoint LINOP of L, LSTAR, under the condition that L is a linear
%   differential operator, i.e. 
%
%       L*u = a_k*diff(u,k) + ... + a_0*u.
%
%   If L represents a system of differential equations then the differential
%   order in each variable must be the same. Integral operators are not
%   supported. 
%
%   The outputs are a LINOP LSTAR that represents the formal adjoint and a
%   function handle OP that can be used to construct a CHEBOP.
%
% See also ADJOINT and ADJOINTBCS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for 2 inputs:
if ( nargin < 2 )
    error('CHEBFUN:LINOP:adjointFormal:inputs', ...
        'ADJOINT requires two inputs.')
end

% Initialize Lstar and op to have correct dimensions:
[nrows, ncols] = size(L);
op = mat2cell(zeros(ncols, nrows), ones(1, ncols), ones(1, nrows));
Lstar = linop(op);

% If L.blocks{ii,jj} is a chebfun convert to operator block:
for ii = 1:nrows
    for jj = 1:ncols
        if ( isa(L.blocks{ii,jj}, 'chebfun') )
            b = L.blocks{ii,jj};
            L.blocks{ii,jj} = operatorBlock.mult(b, b.domain);
        end
    end
end

% Get the domain and the value of the highest derivative:
dom = L.domain;
Lstar.domain = dom;

% Loop through blocks:
op = ';';
argstr = '@(x';
for ii = 1:ncols
    for jj = 1:nrows
        
        % Get local block:
        B = L.blocks{jj,ii};

        % Get the coefficients:
        coeffs = toCoeff(B, pref);
        dor = length(coeffs)-1;

        % Compute the coefficients of the adjoint:
        adjCoeffs = 0*coeffs;
        for k = 0:dor
            for l = 0:k
                adjCoeffs(dor+1-l) = adjCoeffs{dor+1-l} + ...
                    (-1)^k*nchoosek(k,l)*conj(diff(coeffs{dor+1-k}, k-l));
            end
        end

        % Update the blocks of Lstar using adjoint coefficients:
        Lstar.blocks{ii,jj} = 0;
        M = @(f) operatorBlock.mult(f, dom);
        D = @(k) operatorBlock.diff(dom, k);
        for k = 0:dor
            varname = ['a', int2str(ii), int2str(jj), '_', int2str(k)];
            eval([varname, '= adjCoeffs{dor+1-k};']);
            Lstar.blocks{ii,jj} = Lstar.blocks{ii,jj} + ...
				  M(adjCoeffs{dor+1-k}) * D(k);
        end

        % Update argstr:
        if nrows == 1
            vstr = 'v';
        else
            vstr = ['v', int2str(jj)];
        end
        if ( ii == 1 )
            argstr = [argstr, ',', vstr];
        end

        % Set Bop:
        Bop = [];
        astr = ['a', int2str(ii), int2str(jj)];
        for k = 0:dor
            
            acur = [astr, '_', int2str(k)];
            
            % Check acur for 0, 1, or -1:
            acheb = eval(acur);
            if ( norm(acheb) == 0 && dor == 0 )
                acur = '+0*';
            elseif ( norm(acheb) == 0 && dor > 0 )
                acur = [];
            elseif ( length(acheb) == 1 && acheb(0) == 1 )
                acur = '+';
            elseif ( length(acheb) == 1 && acheb(0) == -1 )
                acur = '-';
            else 
                acur = ['+', acur, '.*'];
            end
            
            % Use acur to update op string:
            if ( ~isempty(acur) )
                if ( k == 0 )
                    Bop = [acur, vstr, Bop];
                elseif ( k == 1 )
                    Bop = [acur, 'diff(', vstr, ')', Bop];
                else
                    Bop = [acur, 'diff(', vstr, ',', int2str(k), ')', Bop];
                end
            end
        end
 
        % Update op:
        if ( ~isempty(Bop) && strcmp(op(end), ';') )
            if ( strcmp(Bop(1), '+') )
                Bop = [' ', Bop(2:end)];
            elseif ( strcmp(Bop(1), '-') )
                Bop = [' ', Bop];
            end
        end
        op = [op, Bop];

    end

    % Update op:
    if ( ii < ncols )
        op = [op, ';'];
    end

end 

% Set op:
if ( ncols > 1 )
    op = [argstr, ') [', op(2:end), ' ];'];
else
    op = [argstr, ') ', op(2:end), ';'];
end
eval(['op = ', op]);

end