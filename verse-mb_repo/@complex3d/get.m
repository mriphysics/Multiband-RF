function res = get(params,field)

eval(strcat('res=params.',field,';'));

