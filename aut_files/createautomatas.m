% we need to define the components models based on their DES definitions, this will require to capture all the parts modesl
% such as: Lines, generatoes and loads 
%% step 1 : call function from matpower and define the power system 

% compile the .cpp files for Generator and parallel operation 

mex -lfaudes -L"D:\IEEE_300_Bus\libFAUDES" Automaton.cpp

% prepare automatas for all the lines, loads and generatetors by caling the
% compiled function Automaton.cpp , for each component 



% get the lines, generators and loads indices connected to each bus and do
% parallel operation to these component, bus by bus.  


%% step 2: load case info (mpc)

% mpc = loadcase('case6ww');

mpc = case300_002_ps;


%% step 3: convert the power system buses into Automata using LibFAUDES 

% through mpc.branch get all the transmission lines IDs: 
 
         branchesnumber = size(mpc.branch,1);  % return the number of rows in the matrix mpc.branch (branch map) 
         
% we'll use simpler DES model: Normal and trip states only 
         
        for j = 1:branchesnumber
           
            if j < 10
            x = ['L' '0' '0' num2str(j) 'N', 'L' '0' '0' num2str(j) 'T','line' '0' '0' num2str(j)];
            % then pass the events in a seperate variable 
            y = ['k' '0' '0' num2str(j), 'u' '0' '0' num2str(j), 'h' '0' '0' num2str(j)];
            % calling this function after compiling will generate .gen object file in
            % the same directory 
            Automaton(x,y)  
            else 
                 if j >= 10 && j <= 99

                        x = ['L' '0' num2str(j) 'N', 'L' '0' num2str(j) 'T','line' '0' num2str(j)];
                        % then pass the events in a seperate variable 
                        y = ['k' '0' num2str(j), 'u' '0' num2str(j), 'h' '0' num2str(j)];
                        % calling this function after compiling will generate .gen object file in
                        % the same directory 
                        Automaton(x,y)  

                 else

                        x = ['L' num2str(j) 'N', 'L' num2str(j) 'T','line' num2str(j)];
                        % then pass the events in a seperate variable 
                        y = ['k' num2str(j), 'u' num2str(j), 'h' num2str(j)];
                        % calling this function after compiling will generate .gen object file in
                        % the same directory 
                        Automaton(x,y)  

                 end
            end 
        end 

        % through mpc.bus get all the generators and loads IDs:
% first construt the loads Automata 
Loadbuses = find(mpc.bus(:,3)> 1);         
Loadbusescounter = size(mpc.shunt,1);
       for t =1:Loadbusescounter

            if t < 10
            x = ['D' '0' '0' num2str(t) 'N', 'D' '0' '0' num2str(t) 'T','Load' '0' '0' num2str(t)];
            % then pass the events in a seperate variable 
            y = ['e' '0' '0' num2str(t), 'f' '0' '0' num2str(t), 'g' '0' '0' num2str(t)];
            % calling this function after compiling will generate .gen object file in
            % the same directory 
            Automaton(x,y)  
            else 
                if t >= 10 && t <= 99

                    x = ['D' '0' num2str(t) 'N', 'D' '0' num2str(t) 'T','Load' '0' num2str(t)];
                    % then pass the events in a seperate variable 
                    y = ['e' '0' num2str(t), 'f' '0' num2str(t), 'g' '0' num2str(t)];
                    % calling this function after compiling will generate .gen object file in
                    % the same directory 
                    Automaton(x,y)  

                else

                    x = ['D' num2str(t) 'N', 'D' num2str(t) 'T','Load' num2str(t)];
                    % then pass the events in a seperate variable 
                    y = ['e' num2str(t), 'f' num2str(t), 'g' num2str(t)];
                    % calling this function after compiling will generate .gen object file in
                    % the same directory 
                    Automaton(x,y) 

                end
            end
           
       end 

% Second: construt the Generators Automata

Genbuses = find(mpc.bus(:,2)== 2);         
Genbusescounter = size(mpc.gen,1);  

      for i =1:(Genbusescounter) % reserve bus 1 to inidcate refrence (slack) bus
           
        if i < 10
            x = ['G' '0' '0' num2str(i) 'N', 'G' '0' '0' num2str(i) 'T','Genr' '0' '0' num2str(i)];
            % then pass the events in a seperate variable 
            y = ['a' '0' '0' num2str(i), 'b' '0' '0' num2str(i), 'c' '0' '0' num2str(i)];
            % calling this function after compiling will generate .gen object file in
            % the same directory 
            Automaton(x,y)  
        else 

            if i >= 10 && i <= 99

                x = ['G' '0' num2str(i) 'N', 'G' '0' num2str(i) 'T','Genr' '0' num2str(i)];
                % then pass the events in a seperate variable 
                y = ['a' '0' num2str(i), 'b' '0' num2str(i), 'c' '0' num2str(i)];
                % calling this function after compiling will generate .gen object file in
                % the same directory 
                Automaton(x,y)  

            else 
                x = ['G' num2str(i) 'N', 'G' num2str(i) 'T','Genr' num2str(i)];
                % then pass the events in a seperate variable 
                y = ['a' num2str(i), 'b' num2str(i), 'c' num2str(i)];
                % calling this function after compiling will generate .gen object file in
                % the same directory 
                Automaton(x,y)  

            end



            end   
      end
      
           
      
      


