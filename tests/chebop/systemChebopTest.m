% cheboppref('display','iter')
% cheboppref('plotting','on')
N = chebop(@(x,u,v) [diff(u,2) + sin(v); cos(u) + diff(v,2)]);
N.lbc = @(u,v) [u-2; v-1]; N.rbc = @(u,v) [u - 2, v + 1];
uv = N\0