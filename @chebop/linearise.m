function [L, res, isLinear] = linearise(N, x, u, flag)
isLinear = true(1, 4);
numVars = nargin(N.op) - 1;

if ( nargin < 2 )
    x = chebfun(@(x) x, N.domain);
end
if ( nargin < 3 )
    % Initialise a zero ADCHEBFUN:
    zeroFun = chebfun(0, N.domain);
    u = cell(numVars, 1);
    for k = 1:numVars
        u{k} = zeroFun;
    end
    %                 u0 = chebmatrix(u0);
end

if ( isa(u, 'chebmatrix') )
    u = u.blocks;
end
if ( isa(u, 'chebfun') )
    u = {u};
end

for k = 1:numVars
    u{k} = seed(adchebfun(u{k}), k, numVars);
end


%%
% Evaluate the operators to get a linearisation:

% Experiment with using chebmatrices:
%             w = N.op(x, u{:});
%             if isa(w, 'adchebfun')
%                 L = w.jacobian;
%                 res = w.func;
%                 isLinear(1) = all(w.isConstant);
%             else
%                 L = cellfun(@(w) get(w, 'jacobian'), w.blocks, 'uniformoutput', false);
%                 L = vertcat(L{:});
%                 res = cellfun(@(w) get(w, 'func'), w.blocks, 'uniformoutput', false);
%                 res = chebmatrix(res);
%                 isLinear(1) = all(cellfun(@(w) get(w, 'isConstant'), w.blocks));
%             end

w = N.op(x, u{:});
L = linop(vertcat(w.jacobian));
res = vertcat(w.func);
isLinear(1) = all(w.isConstant);

BC = linopConstraint();
%%
% Add BCs
if ( ~isempty(N.lbc) )
    lbcU = N.lbc(u{:});
    for k = 1:numel(lbcU)
        lbcU(k) = feval(lbcU(k), N.domain(1));
        if ( nargin == 4 ), lbcU(k).func = -lbcU(k).func; end
        BC = append(BC, lbcU(k).jacobian, -lbcU(k).func);
        %                     L = bc(L, lbcU(k).jacobian, -lbcU(k).func); %#ok<CPROP>
    end
    isLinear(2) = all([lbcU.isConstant]);
end

if ( ~isempty(N.rbc) )
    rbcU = N.rbc(u{:});
    for k = 1:numel(rbcU)
        rbcU(k) = feval(rbcU(k), N.domain(end));
        if ( nargin == 4 ), rbcU(k).func = -rbcU(k).func; end
        BC = append(BC, rbcU(k).jacobian, -rbcU(k).func);
        %                     L = bc(L, rbcU(k).jacobian, -rbcU(k).func); %#ok<CPROP>
    end
    isLinear(3) = all([rbcU.isConstant]);
end

if ( ~isempty(N.bc) )
    % Experimant with using chebmatrix to store adchebfuns:
    %                 bcU = N.bc(x, u{:});
    %                 con = cellfun(@(w) get(w, 'jacobian'), bcU.blocks, 'uniformoutput', false);
    %                 con = vertcat(con{:});
    %                 val = cellfun(@(w) get(w, 'func'), bcU.blocks, 'uniformoutput', false);
    %                 val = chebmatrix(val);
    %                 if ( nargin == 4 ), val = -val; end
    %                 L = bc(L, con, -val); %#ok<CPROP>
    
    bcU = N.bc(x, u{:});
    vals = cat(1, bcU.func);
    if ( nargin == 4 ), vals = -vals; end
    for k = 1:numel(bcU)
        BC = append(BC, bcU(k).jacobian, -vals(k));
        %                     L = bc(L, bcU(k).jacobian, -vals(k)); %#ok<CPROP>
    end
    isLinear(4) = all([bcU.isConstant]);
end
L.constraint = BC;

end