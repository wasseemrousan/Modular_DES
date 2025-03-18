function  [TargetNodes,TargetNodesChar] = AutomatonIterator_Matlab_chars(G,s)

% this function has the same job as AutomatonIterator_Matlab_numbers except
% it is for smaller supervisors or grapghs that don't have state numbers in
% the event list , bu the state name directly 

% if s or source node is given in chars then convert it to numbers 



alltargetnodes = distances(G.G,s); % find distances of all nodes from source node s % s here is in chars not numbers 

%% A part for empty result need to be added




        
         ntr = 0; %start counter for desired array of target nodes
        
        %% then find TargetNodes from alltargetnodes with all dsitances 
        
        for y =1:size(alltargetnodes,2) 
        
            if alltargetnodes(y) == 1 % if distance is 1 then there is direct connection and we need this as a target node
        
        % then target node # y is needed 
        
             ntr = ntr +1; 
             TargetNodes(ntr)= y; % this is in the case of all the nodes names is defined in G
             %  TargetNodes is given in chars diretly 
        
            end
        
        
        end 
        
        
            TargetNodesChar = strings(1,size(TargetNodes,2)); % empty string of same size as TargetNodes
            for f = 1:size(TargetNodes,2)
                
               TargetNodesChar(f) =  G.G.Nodes{TargetNodes(f),1};
            
            end




end 