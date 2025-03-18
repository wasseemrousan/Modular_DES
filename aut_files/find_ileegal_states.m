

Busid = {};
for r =1:300
if r <= 9 
    Busid{r} = ['Bus00' num2str(r) 'e.gen'];
end

if r <= 99 && r > 9 
    Busid{r} = ['Bus0' num2str(r) 'e.gen'];
end


if r >= 100
    Busid{r} = ['Bus' num2str(r) 'e.gen'];
end

end

illegalstates = {}; 
initstates= {};
for i=1:300
    
    ille = readlines(Busid{i});
    indxarray = strfind(ille,'<InitStates>');
    indxinit = find(~cellfun(@isempty,indxarray));
    initstates{i}= ille((indxinit+1),1);
    initstates{i} = strtrim(initstates{i});
    illegalstates{i} = ille(indxinit+1);
    illegalstates{i} = strtrim(illegalstates{i});
    illegalstates{i} = replace(illegalstates{i},'N','T');

end