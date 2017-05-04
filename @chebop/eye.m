function I = eye(A)
%EYE   Create an identity CHEBOP
%   I = eye(A) creates an identity CHEBOP with the same domain and
%   size(i.e., number of variables) as the CHEBOP A.

if ( isempty(A) )
    I = chebop(@(u) u);
    return
else
    I = chebop(A.domain);
end

funArgs = getFunArgs(A);
if ( isempty(funArgs) )
    funArgs = 'u'; % Default to u?
end 

idx = strfind(funArgs, ',');
if ( isempty(idx) )
    opStr = funArgs;
elseif ( numel(idx) == 1 )
    args = funArgs;
    opStr = args(idx+1:end);
else
    args = funArgs;
    opStr = [args(idx(end)+1:end) ']'];
    for k = numel(idx):-1:2
        opStr = [args(idx(k-1)+1:idx(k)-1) ';' opStr];
    end
    opStr = ['[' opStr ];
end

I.op = eval(['@(', funArgs, ') ', opStr]); % Create new anon. func

end