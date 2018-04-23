function [root, ea] = newtonraphson (func, xi, es)

if nargin < 2
    error('Function requires 2 inputs')
end 
root = 0;
ea = 100;
while ea > es
    xr = root;
    fp = diff(func);
    root = xi - (f(xi)./fp(xi));
    fprintf('The root for this function is %f\n', root)
    ea = (abs((xr-root)/(root)))*100;
end
end 
