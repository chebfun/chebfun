function C = mtimes(A, B)
%*   Composition of CHEBMATRICES.
%
% See also MPOWER.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Scalar operand is scalar*identity. But it's faster to interpret it as a
% special case. 

if ( isnumeric(A) && (numel(A)==1) )  
    C = scalarTimes(B, A);
    
elseif ( isnumeric(B) && (numel(B)==1) )
    C = scalarTimes(A, B);
    
elseif ( isnumeric(B) )
    C = matrixTimes(A, B);

else
    % Ensure both A and B are CHEBMATRICES so that we can compare blocks:
    if ( ~isa(A, 'chebmatrix') )
        A = chebmatrix(A);
    elseif ( ~isa(B, 'chebmatrix') )
        B = chebmatrix(B);
    end
    
    % Setup before we can start composition:
    [m, n] = size(A);
    Adata = A.blocks;   % Needed to avoid subsref call later.
    Bdata = B.blocks;
    [nB, p] = size(B);
    if ( n ~= nB )
        error('CHEBFUN:CHEBMATRIX:mtimes:dims', ...
            'Operand inner dimensions must agree.')
    end
    
    C = cell(m, p);    
    % Do the composition:
    for i = 1:m
        for j = 1:p
            % It's tricky to just start a sum with "zero", because we don't know
            % if it's added to an operator or a functional.
            u = Adata{i,1} * Bdata{1,j};
            for k = 2:n
                u = u + Adata{i,k} * Bdata{k,j};
            end
            C{i, j} = u;
        end
    end
    % Convert the resulting cell C to a chebmatrix to be returned.
    C = chebmatrix(C);
    
end
    
end

function C = scalarTimes(A, z)
%SCALARTIMES   Multiplying blocks in a CHEBMATRIX with a scalar.
[m, n] = size(A);
C = cell(m, n);
for i = 1:m
    for j = 1:n
        C{i,j} = z * A.blocks{i,j};
    end
end
C = chebmatrix(C);
end

function C = matrixTimes(A, B)
%MATRIXTIMES   Multiply CHEBMATRIX * matrix (or vector).

% All the row-entries of A must be of the same type, so we HORZCAT. This is much
% much much faster than converting a matrix to a CHEBMATRIX and looping over all
% the block entries.
m = size(A,1); 
C = cell(m, size(B,2));
for k = 1:m
    tmp = horzcat(A.blocks{k,:}); % Horzcat
    if ( isa(tmp, 'chebfun') )    % Use array-valued if possible
        tmp = quasi2cheb(tmp);    % Since we're adding, all breaks
    end                           % will eventually be the same. 
    if ( isvector(B) )
        C(k) = tmp*B;             % Avoid num2cell call.
    else
        C(k,:) = num2cell(tmp*B);
    end
end
% Convert cell to CHEBMATRIX for ouput:
C = chebmatrix(C);  
end
