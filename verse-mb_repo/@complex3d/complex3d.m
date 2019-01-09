function res = complex3d(in)

res.adjoint = 0;
res.PP =[];

if nargin == 0
    res.data = zeros(0,3);
    res = class(res,'complex3d');
elseif isa(in,'complex3d')
    res = in;
else
    res.data = in;
    res = class(res,'complex3d');
end
