function res = mtimes(a,b)

if isa(a,'complex3d')
    if isa(b,'complex3d')
        data = a.data*b.data;
    else
        data = a.data*b;
    end
else
    data = a*b.data;
end

res = complex3d(data);


