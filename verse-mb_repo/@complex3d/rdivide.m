function res = rdivide(a,b)

if isa(a,'complex3d')
    if isa(b,'complex3d')
        data = a.data./b.data;
    else
        data = a.data./repmat(b,[1,3]);
    end
else
    data = a./b.data;
end

res = complex3d(data);


