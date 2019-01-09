function params = set(params,field,val)


p = struct(params);

if isfield(p,field)
	eval(strcat('params.',field,'= val;'));
else
	error('no Such field');
end

