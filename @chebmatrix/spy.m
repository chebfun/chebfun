function spy(A)
%SPY    Visualize a chebmatrix.
%   If A is a chebmatrix, SPY(A) creates a picture of the nonzero pattern of the
%   default discretization of A.
%   
%   SPY(A,DIM) uses the dimension vector DIM to create the picture.
%
%   See also CHEBMATRIX, CHEBMATRIX.MATRIX, CHEBOPPREF.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Select a reasonable default discretization.
dom = A.domain;
if length(dom)==2   % no breakpoints
    dim = 10;
else
    dim = repmat(6,[1 length(dom)-1]);
end

% Discretize and do a regular spy plot.
data = matrix(A,dim,dom);
spy(data)
s =  sprintf('%i,',dim);    % list of sizes
s = [ 'discretization = [', s(1:end-1), ']' ];
xlabel(s)

% Reverse engineer to see where the block boundaries are.
[m,n] = blockSizes(A);
m(isinf(m)) = sum(dim);   % use discretization sizes
n(isinf(n)) = sum(dim);
csrow =cumsum(m(:,1)');
rowdiv = csrow(1:end-1) + 1/2;   % insert boundary after each block
cscol =cumsum(n(1,:));
coldiv = cscol(1:end-1) + 1/2; 

hs = ishold;
hold on
plot([rowdiv;rowdiv],[0;cscol(end)+1],'color',[.6 .6 .6])
plot([0;csrow(end)+1],[coldiv;coldiv],'color',[.6 .6 .6])
if ( ~hs )
    hold off
end

end
