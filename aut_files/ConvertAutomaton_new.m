function  [G,Nodes_Number,Nodes_Names] = ConvertAutomaton_new(G1)

% this is used for supervisor18.gen file size 

%     A = fscanf(fopen('line25.gen'),'%s')
A1 = readlines(G1); % maybe include white paces with %c for better alignement on transistions
% definw states s
clear G1;
             % find number of states 
             maxstringsize = size(A1{6},2); 
             numstates = A1{6}(15:maxstringsize); % how to handlw this if number of states is more than 9 ? same question for number of events ? maybe solution is to make new function and call any of them if needed from main window
             statesnumbers = str2double(numstates);
             % start geeting states names
             r=0;
             k=0;
%              states = zeros(1:statesnumbers);
             states= strings([1,statesnumbers]);
             % Id= strings([1,statesnumbers]);
%              for t=1:statesnumbers
%                  
%                  states{t} = A1(e+8+r:e+224+r); 
%              r = r+226; % maybe this needs to be more dynamic and takes size of chars of states
%              end 
%              names = states;
             names = 1:statesnumbers; %in very large automatas libfaudes start to give states numbers
             % a = zeros(1,size(e,2));

             % start a for loop to get states names in char rather than
             % numbers only 
             e = strfind(A1,'<States>');
             e2= find(~cellfun(@isempty,e));
             e3 = strfind(A1,'</States>'); % find begining of states id index
             e4= find(~cellfun(@isempty,e3)); % find end of the states ids index
             
             %e2 = e2 +
             r = 1;
             % nodenames = strings(1,statesnumbers);
             for u =(e2+1):(e4-1)
                
                TF = isspace(A1{u});
                space_id = find(TF);
                for g = 1:size(space_id,2)
                    if g ==1 
                            nodenames{r} = A1{u}(1:space_id(g)); 
                            e5 = strfind(nodenames{1,r},'#'); % find begining of states id index
                            maxstringsize = strlength(nodenames(r)); 
                            Id{r}= nodenames{r}((e5+1):maxstringsize); % find end of the states ids index
                            Id{r} = str2double(Id{1,r});
                            r = r + 1;
                    else
                        if ((space_id(g)) - space_id(g-1)) ~= 1 
                            nodenames{r} = A1{u}(space_id(g-1):space_id(g)); 
                            e5 = strfind(nodenames{1,r},'#'); % find begining of states id index
                            Id{r}= nodenames{r}((e5+1):maxstringsize); % find end of the states ids index
                            Id{r} = str2double(Id{1,r});
                            r = r + 1;
                        end
                    
                    end
                    clear e5;
                 
                end
             
             end
             clear TF;
             clear space_id;
             clear e;
             clear e2;
             clear e3;
             clear e4;

    % Now define events
    eventsloc = strfind(A1,'<TransRel>');
    events= find(~cellfun(@isempty,eventsloc));
    clear eventsloc;
    maxstringsize_events = size(A1{9},2); 
    evnumloc = A1{9}(16:maxstringsize_events);
    numberofevents = str2double(evnumloc);

    
    
    sourcetstates = strings([1,numberofevents]);
    endstates = strings([1,numberofevents]);
    weights = strings([1,numberofevents]);

    


    for w =1:numberofevents
    TF_events = isspace(A1{events+w});
    space_id__events = find(TF_events);
   
    maxstringsize_state = size(A1{events+w},2); 
    sourcetstates(w) = A1{events+w}(1:space_id__events(1));
    sourcetstates(w) = str2double(sourcetstates(w));

    j=1;   
    for z= 2:size(space_id__events,2)
            if (((space_id__events(z)) - space_id__events(z-1)) ~= 1) && (j==1)
            % define event names
                EdgesName{w} = A1{events+w}(space_id__events(z-1)+1:(space_id__events(z-1)+4));

            j = j+1; 
            end
            
            if (((space_id__events(z)) - space_id__events(z-1)) ~= 1) && (j==2)
                endstates(w) = A1{events+w}(space_id__events(z):maxstringsize_state);
                endstates(w) = str2double(endstates(w));
                j = j+1; 
            end  
      
    end
                     

    end
    clear A1;

   disp('Convert Imported Supervisor Automaton from LIBFaudes into a matlab directed graph')           
    
                                      
    sourcetstatesdouble = str2double(sourcetstates);
    endstatesdouble = str2double(endstates);         
%     weights = [1:6780];
    
    % G = digraph(sourcetstates,endstates,weights,names);
%     G = digraph(sourcetstatesdouble,endstatesdouble);
%     G.Edges.Name = EdgesName';
%     G.Edges.Weight = ones(size(EdgesName,2),1);
%     G = digraph(sourcetstatesdouble,endstatesdouble,ones(size(EdgesName,2),1));
    % Add nodes names here
    EdgeTable = table([sourcetstatesdouble' endstatesdouble'],ones(size(EdgesName,2),1),EdgesName', ...
    'VariableNames',{'EndNodes' 'Weight' 'EdgesName'});
    G = digraph(EdgeTable); 
%     
%     B = convertStringsToChars(nodenames);
% %     finalnames = num2cell(B);
%     G.Nodes.Names = B';
    Id = cell2mat(Id); 
    for o =1:size(Id,2)
        if  Id(o) < o
            Id(o) = o;
        end
    end
    
    
    Nodes_Number = Id';
    Nodes_Names = nodenames'; 
    % define adjacency matrix 
%     A = adjacency(G);  % A alos uses much space
%     B = full(A); This will make memory issues , lets use A for now
%     G.adjacency = A;
%     G.full = B;
end 