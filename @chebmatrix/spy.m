function output = spy(A)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
dom = A.domain;
if length(dom)==2
    dim = 12;
else
    dim = repmat(8,[1 length(dom)-1]);
end
data = matrix(A,dim,dom);

h = cellplot(data);
% CELLPLOT seems to cover up the text representations of double
% values. We give a positive z-value so that they sit on top again.
% And we hide its big ugly box.
for i = 2:length(h)
    if strcmp(get(h(i), 'type'), 'text')
        set(h(i), 'position', [0 0 1]+get(h(i), 'position'))
        set(h(i-1), 'vis', 'off')
    end
end
if ( nargout > 0 )
    output = h;
end

end
