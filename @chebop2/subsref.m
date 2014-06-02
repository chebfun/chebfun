function varargout = subsref(N, ref)
% SUBSREF for CHEBOP2
%
% N.PROP  return the CHEBOP2 property PROP of N. This is equivalent to 
% GET(N, 'PROP'). 
% 
% N(m, n) returns a mn x mn matrix representing a discretization of N. 
% 
% N( F ) return the forward application of N to the chebfun2 F. This is
% equivalent to N * F. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

indx = ref(1).subs;

switch ( ref(1).type )
    case '.'
        if numel(ref) == 1
            % This is a GET call for a property:
            varargout = { get(N,indx) };
        end
    case '()'
        % return discretisation of size indx{1} x indx{2}. 
        if ( length(indx) == 2 )
            n=indx{1}; m=indx{2};
            CC = chebop2.constructDiscretisation(N, chebfun2(0), m+2, n+2);
            sz = size(CC{1,1},1)*size(CC{1,2},2);
            
            % return as a matrix: 
            A = spalloc(sz, sz, m*sz + sz);
            for jj = 1:size(CC,1)
                A = A + kron(CC{jj,1},CC{jj,2});
            end
            
            % If the matrix is not too large, convert to full: 
            if max(size(A))>100
                varargout = { A }; 
            else
                varargout = { full( A ) };
            end
        elseif ( length(indx) == 1 && isa(indx{1},'double') )
            % return a square discretization: 
            ref.type = '()'; 
            ref.subs = {indx{1},indx{1}}; 
            varargout = {subsref(N, ref)};  
        elseif ( length(indx) == 1 && isa(indx{1},'chebfun2') )
            % Forward application of CHEBOP2 to a CHEBFUN2 
            op = N.op; 
            varargout = { op(indx{1}) }; 
        else
            error('CHEBFUN2v:subsref:nonnumeric',...
                'nonnumeric value is not recognised.')
        end
    otherwise
        error('CHEBFUN2v:UnexpectedType',['??? Unexpected index.type of ' index(1).type]);
end