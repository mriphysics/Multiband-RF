function a = subsref(a,s)

if strcmp(s.type,'()')
    if length(s.subs) > 1
        error('Error, susbs is more than 1');
    end
    a.data = a.data(s.subs{1},:);
else
    error('doesnt support any other parentheses than () ');
end


