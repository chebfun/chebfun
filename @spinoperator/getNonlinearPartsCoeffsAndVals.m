 function [Nc, Nv] = getNonlinearPartsCoeffsAndVals(S)
 %GETLINEARANDNONLINEARPARTS   Get the nonlinear parts in coefficient and
 %value space of a 1D, 2D or 3D time-dependent PDE.
 %   [Nc, Nv] = GETNONLINEARPARTSCOEFFSANDVALS(N), where S is a SPINOPERATOR,
 %   outputs two function handles NC and NV which represent the nonlinear parts
 %   in coefficient and value space, of the PDE represented by S.
 
 % Get the variables of the workspace:
 N = S.nonlinearPart;
 func = functions(N);
 wrk = func.workspace{1};
 names = fieldnames(wrk);
 if ( isempty(names) == 0 )
     lengthNames = size(names, 1);
     for k = 1:lengthNames
         eval(sprintf('%s = wrk.(names{k});', names{k}));
     end
 end
 
 % Convert FUNCN to a string STRN:
 strN = func2str(N);
 
 % Get the number of variables NVARS:
 nVars = nargin(N);
 
 % Get the dimension DIM:
 dim = S.dimension;
 
 % For scalar equations in 1D, we support nonlinearities of the form
 % diff(f(u),m) with m>=0:
 if ( dim == 1 && nVars == 1 )
     
     % Get rid of the differentiation part in STRN to get NV=f(u):
     oldString = {'diff', ',\d*)'};
     newString = {'', ')'};
     Nv = regexprep(strN, oldString, newString);
     Nv = eval(Nv);

     % Compute the differentiation order to get NC=diff(u,m):
     diffOrderTwoOrGreater = regexp(strN, ',\d*', 'match');
     diffOrderOne = regexp(strN, 'diff(', 'match');
     if ( isempty(diffOrderTwoOrGreater) == 0 )
         diffOrder = diffOrderTwoOrGreater{1}(2:end);
     elseif ( isempty(diffOrderOne) == 0 )
         diffOrder = '1';
     else
         diffOrder = '0';
     end
     Nc = ['@(u) diff(u,', diffOrder, ')'];
     Nc = str2func(Nc);
     
 % For scalar equations in 2D and 3D, and for systems of equations in 1D, 2D and
 % 3D, we only support nonlinearities of the form f_i(u_1,...,u_n), i.e., with 
 % no differentiation, so Nv=N and Nc=1:
 else
     
     % Nv=N and Nc=1:
     strNv = func2str(N);
     Nc = 1;
     
     % We're going to relabel the variables, e.g., @(u,v) u + v is going
     % to be relabelled @(u) u(1:length(u)/2) + u(length(u)/2+1:end).
     % First, get the names of the variables:
     openParenthesis = strfind(strNv, '(');
     openParenthesis = openParenthesis(1);
     closeParenthesis = strfind(strNv, ')');
     closeParenthesis = closeParenthesis(1);
     variablesNames = strNv(openParenthesis+1:closeParenthesis-1);
     variablesNames = regexp(variablesNames,  ',', 'split');
     
     % Second, relabel the variables:
     strNv = strNv(closeParenthesis+1:end);
     for k = 1:nVars
         idx1 = [num2str((k-1)/nVars), '*', 'length(u)', '+1'];
         idx2 = [num2str(k/nVars), '*', 'length(u)'];
         strNvNew = strrep(strNv, variablesNames{k}, ...
             ['u(', idx1, ':', idx2, ')']);
         strNv = strNvNew;
     end
     Nv = str2func(['@(u)', strNvNew]);
     
 end
 
 end