function a = subsasgn(a,s,b)

if strcmp(s.type,'()')
    if length(s.subs) > 1
        error('Error, susbs is more than 1');
    end
    a.data(s.subs{1},:) = b.data;
else
    error('doesnt support any other parentheses than () ');
end


