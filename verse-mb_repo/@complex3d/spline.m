function a = spline(p,a)

p  = p(:);

a.PP  = {spline(p,a.data(:,1)), spline(p,a.data(:,2)), spline(p,a.data(:,3))};
a.data = [];


