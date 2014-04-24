function varargout = equationPoints(disc)
%EQUATIONPOINTS   Points at which collocation is enforced.

% In COLLOC2, functions are discretized at 2nd kind points but equations are
% enforced at 1st kind points, to avoid duplication at boundaries.

[varargout{1:nargout}] = colloc.points(disc, 1);

end