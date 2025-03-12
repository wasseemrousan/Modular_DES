
load initstates;
load illegalstates;

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



for i=1:300

%% how to read illegal states in automatic way ??

%%
S = readlines('makespec_general_300bus.cpp');
strnum= strlength(illegalstates{1,i});
% find the line number 81 and change number of chars in illegal states
numberofchars= strnum+1;
S{81} =     ['char illegal[' num2str(numberofchars) '];'];
S{83} =     ['for (int i = 0; i <' num2str(numberofchars) '; i++)'];
S{85} =     ['if (i !=' num2str(strnum) ')'];

[fid, msg] = fopen('makespec_general_300bus.cpp', 'w'); 

if fid < 1;
error('could not write output file because "%s"', msg);
end
fwrite(fid, strjoin(S, '\n'));
fclose(fid);

mex -lfaudes -L"D:\IEEE_300_Bus\libFAUDES" makespec_general_300bus.cpp

makespec_general_300bus(Busid{1,i},convertStringsToChars(illegalstates{1,i})); 


% now change the file name of the specifcaition 
if i <= 9 
      newname{i} = ['Bus' '00' num2str(i) 'spec.gen']; 
else 
      if i > 9 && i <= 99
    
          newname{i} = ['Bus' '0' num2str(i) 'spec.gen'];
      else 
          newname{i} = ['Bus' num2str(i) 'spec.gen'];
      end
end

movefile('tmp_specification12.gen',newname{i})


end