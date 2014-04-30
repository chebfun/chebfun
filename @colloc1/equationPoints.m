function varargout = equationPoints(disc)
%EQUATIONPOINTS   Points at which collocation is enforced.

[varargout{1:nargout}] = colloc.points(disc, 1);

end