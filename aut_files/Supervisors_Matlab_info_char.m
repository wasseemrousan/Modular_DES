

for t =10:10

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

   
end