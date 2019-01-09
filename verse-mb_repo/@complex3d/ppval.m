function a = ppval(a, p)

p = p(:);

a.data = [ppval(a.PP{1},p), ppval(a.PP{2},p) , ppval(a.PP{3},p) ];
a.PP = [];

