function varargout = subsref(N,ref)
% SUBSREF for chebop2.

% Copyright 2012 by The University of Oxford and The Chebfun2 Developers.

indx = ref(1).subs;

switch ref(1).type
    case '.'
        if numel(ref) == 1
            % This is a get call to get a property.
            varargout = { get(N,indx) };
        end
    case '()'
        % return discretisation...
        if length(indx) == 2
            n=indx{1}; m=indx{2};
            CC=constructDiscretisation(N,chebfun2(0),m+2,n+2);
            sz = size(CC{1,1},1)*size(CC{1,2},2);
            A = spalloc(sz, sz, m*sz + sz);
            for jj = 1:size(CC,1)
                A = A + kron(CC{jj,1},CC{jj,2});
            end
            if max(size(A))>100
                varargout = {A}; 
            else
                varargout = {full(A)};
            end
            return;
        elseif length(indx) == 1 && isa(indx{1},'chebfun2')
            op = N.op; 
            varargout = { op(indx{1}) }; 
            return
        end
        error('CHEBFUN2v:subsref:nonnumeric',...
            'nonnumeric value is not recognised.')
    otherwise
        error('CHEBFUN2v:UnexpectedType',['??? Unexpected index.type of ' index(1).type]);
end