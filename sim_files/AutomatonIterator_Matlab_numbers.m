function  [TargetNodes, TargetNodesChar] = AutomatonIterator_Matlab_numbers(G,s)

% 504 stands for number of states in the .gen file that we are converting 
% [___] = allpaths(G,s,t,Name,Value) specifies additional options using one or more name-value arguments. 
% You can use any of the output argument combinations in previous syntaxes. For example, 
% you can specify MaxNumPaths and a scalar to limit the number of paths returned

%***************************************************%
%  paths1 = allpaths(G,2,5,'MaxNumPaths',10)
%***************************************************%

%% 
% graph search methods ???? 

% or shall we use dfsearch
% Depth-first graph search ?????????


% or bfsearch  
% Breadth-first graph search  ???? 


%% I think we shall use the funstion distances(G) which is better or distances(G,s,t)

% in this case we need to define state by its number only !!! example
% distances(G,1), which first node as defined in G.nodes ?
% We are only intrestd with the nodes from source node s 
%alltargetnodes = distances(G,1); % find distances of all nodes from source node s
% Wasseem 10/03/2022 ; s shall be a number 
alltargetnodes = distances(G.G,s); % find distances of all nodes from source node s

ntr = 0;%start counter for desired array of target nodes

%% then find TargetNodes from alltargetnodes with all dsitances 

for y =1:size(alltargetnodes,2) 

    if alltargetnodes(y) == 1 % if distance is 1 then there is direct connection and we need this as a target node

% then target node # y is needed 

     ntr = ntr +1; 
     TargetNodes(ntr)= y; % this is in the case of all the nodes names is defined in G
  % !!! we also want the actual value [char values of the target node ] not only the number 
  % how to map node number to 

    end


end 


    TargetNodesChar = strings(1,size(TargetNodes,2)); % empty string of same size as TargetNodes
    cntr = 0; 
    for f = 1:size(G.Nodes_Number,1)
        for d=1:size(TargetNodes,2)
    
            if G.Nodes_Number(f) == TargetNodes(d)
                cntr = cntr +1;
                TargetNodesChar(cntr) = G.Nodes_Names(f);
             
            end
            
        end
    
    end
end 