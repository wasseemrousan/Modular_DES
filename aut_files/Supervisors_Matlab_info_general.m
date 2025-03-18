


for t =268:268
    if t <= 9 
        Supname{t} = ['Bus' '00' num2str(t) 'sup.gen']; 
    else 
          if (t >= 9) && (t <= 99)
              Supname{t} = ['Bus' '0' num2str(t) 'sup.gen'];
          else
              Supname{t} = ['Bus' num2str(t) 'sup.gen'];
          end
    end


    A1 = readlines(Supname{t}); % maybe include white paces with %c for better alignement on transistions
    v = strfind(A1,'<TransRel>');
    cellindx = find(~cellfun(@isempty,v));
    events_table_type = A1{cellindx+1}; 
    events_table_type = strtrim(events_table_type); % reomve white space
    events_table_type =erase(events_table_type," "); % remove spaces between the chars
    L = strlength(events_table_type); 
    if L > 20 % then this is char type event table => use ConvertAutomaton_new_char function
    % now we need to check
    
      if t <= 9 
        Supname{t} = ['Bus' '00' num2str(t) 'sup.gen']; 
        [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new_char(Supname{1,t}); 
        save(['Bus' '00' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');
      else 
          if (t >= 9) && (t <= 99)
              Supname{t} = ['Bus' '0' num2str(t) 'sup.gen'];
              [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new_char(Supname{1,t}); 
              save(['Bus' '0' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');

          else
              Supname{t} = ['Bus' num2str(t) 'sup.gen'];
              [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new_char(Supname{1,t}); 
              save(['Bus' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');

          end
      end
    
    else

      if t <= 9 
        Supname{t} = ['Bus' '00' num2str(t) 'sup.gen']; 
        [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new(Supname{1,t}); 
        save(['Bus' '00' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');
      else 
          if (t >= 9) && (t <= 99)
              Supname{t} = ['Bus' '0' num2str(t) 'sup.gen'];
              [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new(Supname{1,t}); 
              save(['Bus' '0' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');

          else
              Supname{t} = ['Bus' num2str(t) 'sup.gen'];
              [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new(Supname{1,t}); 
              save(['Bus' , num2str(t), 'supData.mat'], 'G','Nodes_Number','Nodes_Names', '-v7.3');

          end
      end

    end
end