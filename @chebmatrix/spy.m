function output = spy(A)

dom = A.domain;
if length(dom)==2
    dim = 12;
else
    dim = repmat(8,[1 length(dom)-1]);
end
d = A.discretizer(A,dim,dom);
data = matrix(d);
%data = discretizeBlocks(A, 10);
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
