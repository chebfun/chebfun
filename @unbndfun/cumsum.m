function g = cumsum(g)
% CUMSUM	Indefinite integral

% Chain rule:
rescale = onefun.constructor(@(x) g.mapping.forder(x)); % [TODO]: singfun is needed here!
g.onefun = g.onefun.*rescale;
g.onefun = cumsum(g.onefun);

end
