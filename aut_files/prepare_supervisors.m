clc
clear all

load illegalstates.mat;
load initstates.mat;
% convert illegal states from string to chars
for t =1:300 
    illegalstates{1,t} = convertStringsToChars(illegalstates{1,t}); 
end

Busid = {};
specid = {}; 
for r =1:300
if r <= 9 
    Busid{r} = ['Bus00' num2str(r) 'e.gen'];
    specid{r} = ['Bus00' num2str(r) 'spec.gen']; 
end

if r <= 99 && r > 9 
    Busid{r} = ['Bus0' num2str(r) 'e.gen'];
   specid{r} = ['Bus0' num2str(r) 'spec.gen']; 

end


if r >= 100
    Busid{r} = ['Bus' num2str(r) 'e.gen'];
   specid{r} = ['Bus' num2str(r) 'spec.gen']; 

end

end




for i =1:300


S = readlines('syn_6_compsynth.cpp');
% S{134} = regexprep(S{134}, 'contevents.Insert("k19");', 'contevents.Insert("k20");');    %file says -1 change it to 31419
% insert contralble evnts here 
% 1st check for lines ids, gen ids, loads ids 

% Find Gen the substring in.
y = strfind(illegalstates{1,i},'G');
GenCont = {};
for b =1:size(y,2)
GenCont{b} = ['contevents.Insert(','"b',illegalstates{1,i}(y(b)+1),illegalstates{1,i}(y(b)+2),illegalstates{1,i}(y(b)+3),'");']; 
end

% Find Load the substring in.
a = strfind(illegalstates{1,i},'D');
LodCont = {};
for b =1:size(a,2)
LodCont{b} = ['contevents.Insert(','"f',illegalstates{1,i}(a(b)+1),illegalstates{1,i}(a(b)+2),illegalstates{1,i}(a(b)+3),'");']; 
end

% Find Lines the substring in.
f = strfind(illegalstates{1,i},'L');
lineCont = {};
for b =1:size(f,2)
lineCont{b} = ['contevents.Insert(','"k', illegalstates{1,i}(f(b)+1),illegalstates{1,i}(f(b)+2),illegalstates{1,i}(f(b)+3) , '");']; 
end
S{125} = cell2mat([ GenCont, LodCont, lineCont]);

%S{131} = ['contevents.Insert("k20");contevents.Insert("k21");contevents.Insert("k22");contevents.Insert("k23");'];

[fid, msg] = fopen('sup_synth_new.cpp', 'w'); 

if fid < 1;
error('could not write output file because "%s"', msg);
end
fwrite(fid, strjoin(S, '\n'));
fclose(fid);
mex -lfaudes -L"D:\IEEE_300_Bus\libFAUDES" sup_synth_new.cpp
sup_synth_new(Busid{1,i},specid{1,i});




if i <= 9 
      newname{i} = ['Bus' '00' num2str(i) 'sup.gen']; 
else 
      if i > 9 && i <= 99
    
          newname{i} = ['Bus' '0' num2str(i) 'sup.gen'];
      else 
          newname{i} = ['Bus' num2str(i) 'sup.gen'];
      end
end

movefile('tmp_supervisor12.gen',newname{i})


end
