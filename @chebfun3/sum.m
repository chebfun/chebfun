function ff = sum(f, dim)
%SUM   Definite Integration of a CHEBFUN3.
%
%   G = sum(F,DIM) where DIM is 1, 2 or 3 integrates only over X or Y or Z 
%   respectively, a CHEBFUN2 in the remaining variables.
%
%   G = sum(F) is the same as sum(F,1)
%
%   See also chebfun3/sum2, chebfun3/sum3, chebfun3/cumsum, 
%   chebfun3/cumsum2 and chebfun3/cumsum3.

% Empty check: 
if ( isempty(f) ) 
    ff = []; 
    return; 
end

% Default to y direction: 
if ( nargin == 1 )
    dim = 1;
end

% Get the low rank representation for f. 
% cols = f.cols;
% rows = f.rows;
% tubes = f.tubes;

%[cols, D, rows] = cdr(f);
%dom = f.domain;

% For now, suppose fibers = rows, i.e., f.fiberDim = 1.
if ( dim == 1 )
%if ( dim == f.fiberDim )
    % Integrate over x: 
    core = squeeze(chebfun3.txm(f.core, sum(f.cols), dim));
    %ff = chebfun2(f.rows*core*(f.tubes)'); % order is wrong.
    ff = chebfun2(f.tubes*core'*f.rows');
%     ff = chebfun2();
%     [U,S,V] = svd(core); % TODO: replace svd by aca.
%     ff.pivotValues = diag(1./S);
%     ff.rows = f.rows*U;
%     ff.cols = f.tubes*V;
%     %ff = chebfun2(ff.rows.*ff.pivotValues.*ff.cols)
%    ff.domain = [dom(3) dom(4) dom(5)  dom(6)];
%    ff = chebfun2(ff);
    if ( isa(ff, 'chebfun2') ) 
        ff = simplify(ff); 
%     else
%         % f = double 
%         ff = chebfun2(f, dom(3:6)); % What does this mean ??? 
    end
elseif ( dim == 2 )
%elseif ( abs(dim - f.fiberDim) == 1 )
    % Integrate over y: 
    core = squeeze(chebfun3.txm(f.core, sum(f.rows), dim));
    ff = chebfun2(f.tubes*core'*f.cols');
%     ff = chebfun2();
%     [U,S,V] = svd(core);
%     ff.pivotValues = diag(1./S);
%     ff.rows = f.cols*U;
%     ff.cols = f.tubes*V;
%     ff.domain = [dom(1) dom(2) dom(5)  dom(6)];
%     ff = chebfun2(ff);
    if ( isa(ff, 'chebfun2') ) 
        ff = simplify(ff);
    end
    elseif ( dim == 3 )
    % Integrate over z:  
    core = squeeze(chebfun3.txm(f.core, sum(f.tubes), dim));
    ff = chebfun2(f.rows*core'*f.cols');
%     ff = chebfun2();
%     [U,S,V] = svd(core);
%     ff.pivotValues = diag(1./S);
%     ff.rows = f.cols*U;
%     ff.cols = f.rows*V;
%     ff.domain = [dom(1) dom(2) dom(3)  dom(4)];
%     ff = chebfun2(ff);
    if ( isa(ff, 'chebfun2') )
        ff = simplify(ff); 
    end
else 
    error('CHEBFUN:CHEBFUN3:sum:unknown', ...
          'Undefined function ''sum'' for that dimension');
end

end
