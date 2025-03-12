function [TargetNodesChar] = SupervisorIterator(supervisor,startstate,supervisorname)


% A condition needs to be added if we have to use Automaton_Iterator40 or
% Automaton_Iterator60 for bus 2 and 5 

 

if ~isempty(strfind(supervisor.G.Edges.EndNodes(1),'|'))  % if the char elemnt '|' exist in the edges table to see what type of itertator to use
   
             % startstate_index = strfind(supervisor.G.Nodes,startstate);
             startstate_index = find(~cellfun(@isempty,(strfind(supervisor.G.Nodes{:,1},startstate))'));

                State_index = 0; 
                % for j =1:size(startstate_index,1)
                % TF = isempty(startstate_index{j,1});
                %     if TF == 0
                %         State_index = supervisor.Nodes_Number(j);
                %     end
                % end               
                State_index = supervisor.Nodes_Number(startstate_index);
                disp('it is desired to iterate supervisor for bus number :   ')
                disp(supervisorname)
                
                if State_index == 0
                TargetNodesChar = []; 
                else
                [TargetNodes,TargetNodesChar] = AutomatonIterator_Matlab_chars(supervisor,State_index);
                end
                                                        

else

                startstate_index = strfind(supervisor.Nodes_Names,startstate);
                State_index = 0;
                for j =1:size(startstate_index,1)
                TF = isempty(startstate_index{j,1});
                    if TF == 0
                        State_index = supervisor.Nodes_Number(j);
                    end
                end 
              disp('it is desired to iterate supervisor for bus number :   ')
              disp(supervisorname)
              
              if State_index == 0
              
                  TargetNodesChar = []; 
              else
              
                  [TargetNodes, TargetNodesChar] =  AutomatonIterator_Matlab_numbers(supervisor,State_index);
                            
              end            
end