clc
% clear mex
bdclose all
clear all
%% this is the part to:
% 1. build the automatos of the power sysytem, differnt m file than the
% controller 
% 2. send the currrent state string to the controller 
% 3. receive the control actions from the controller
% also make sure that DCSimsep library is added to bath

%% get constants that help us to find the data
C = psconstants; % tells me where to find my data

%% set some options
opt = psoptions;
opt.verbose = false; % set this to false if you don't want stuff on the command line
% Stopping criterion: (set to zero to simulate a complete cascade)
opt.sim.stop_threshold = 0.5; % the fraction of nodes, at which to declare a major separation
% opt.sim.fast_ramp_mins = 1;
opt.sim.simple_redispatch = false;
opt.sim.stop_on_sep = true; 
%% Prepare and run the simulation for the Polish grid

ps = case300_002_ps;

psL = ps; % use psL for load variations purpose 

fprintf('----------------------------------------------------------\n');
disp('loading the data');
tic

ps = case300_002_ps;

toc
fprintf('----------------------------------------------------------\n');
tic
ps = updateps(ps);
ps = rebalance(ps);
ps = dcpf(ps);
toc
fprintf('----------------------------------------------------------\n');
m = size(ps.branch,1);
pre_contingency_flows = ps.branch(:,C.br.Pf);
phase_angles_degrees = ps.bus(:,C.bu.Vang);
Pd_total = sum(ps.shunt(:,C.sh.P));
% Set lower gen limits to zero
ps.gen(:,C.ge.Pmin) = 0;
%%
opt.verbose = false;

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

NL = 10; % this is acounter for the loads variations selection loop

%% define some nesecarry variables for the LLP DES controller 

for b = 1:size(ps.branch(:,1),1)
      
    if b < 10
        currentstatesL{1, b} = ['L' '00' num2str(b) 'N'];
    else 
        if (b >= 10) &&  (b <= 99)
        currentstatesL{1, b} = ['L0' num2str(b) 'N'];
        else
        currentstatesL{1, b} = ['L' num2str(b) 'N'];

        end
    end
       
end

% generators states 
for a = 1:size(ps.gen(:,1),1)
    if a < 10
        currentstatesG{1, a} = ['G' '00' num2str(a) 'N'];
    else 
                if (a >= 10) &&  (a <= 99)

                    currentstatesG{1, a} = ['G0' num2str(a) 'N'];
                else

                    currentstatesG{1, a} = ['G' num2str(a) 'N'];
                end
    end
end

% Load states 
for c = 1:size(ps.shunt(:,1),1)
    if c < 10
        currentstatesD{1, c} = ['D' '00' num2str(c) 'N'];
    else 
        if (c >= 10) &&  (c <= 99)
                    currentstatesD{1, c} = ['D0' num2str(c) 'N'];
        else
                    currentstatesD{1, c} = ['D' num2str(c) 'N'];

        end
    end
end

currentstates_LLP = [currentstatesL currentstatesG currentstatesD];
load initstates;
load illegalstates;
currentstates = initstates;
for i=1:300
currentstates{1,i} = convertStringsToChars(currentstates{1,i});
end
global CrticalLineIndex; 
global alcritical;

MU_Loads = [ps.shunt(:,2)]'; % determine the means of the loads 

global f_new;
global alIndex; 
global alvalue;
global pattern;
global sigmaFl;
global uFt;



%% +++++++++++++++++++++++++++++++10++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% generate 1000 samples of loads and use them for all cases, we will need to

MU_Loads = [ps.shunt(:,2)]'; % determine the means of the loads 
for i = 1:1000
p{i,1} = randperm(411,2);

    for k = 1:size(ps.shunt,1)
        if MU_Loads(k) > 0
         psloads{1,i}(k) = normrnd(MU_Loads(k),(0.15*MU_Loads(k)));
        else 
         psloads{1,i}(k) = -normrnd(abs(MU_Loads(k)),(0.15*abs(MU_Loads(k))));
        end
    end
p{i,2} = psloads{1,i};

end


%main loop for monte carlo (MC) simulation
for i = 1:100 
    
   % the current state string needs to be modifed to match the branch outage continegency 
   currentstatefortest = currentstates;
   BOPAIRS2(i,1) = p{i,1}(1); % first line intentianally tripped 
   BOPAIRS2(i,2) = p{i,1}(2); % first line intentianally tripped
   
%% Prepare some daat for LLP method   
  % the current state string needs to be modifed to match the branch outage continegency 
   BOPAIRS2(i,1); % first line intentianally tripped 
   BOPAIRS2(i,2); % first line intentianally tripped

     
%% Start a for loop for modular supervisors to lable initial state with the contingency,

   if BOPAIRS2(i,1)<10 
   idx1 = strfind(currentstatefortest,['L','00',num2str(BOPAIRS2(i,1)),'N']);
   % because states are repeated in supversiors automaton
        for z = 1:300
        TF = isempty(idx1{1,z});
        if TF == 0 
        currentstatefortest{1,z}(idx1{1,z}+4) ='T';
        end
        end   
     
          
   else
       if (BOPAIRS2(i,1)>=10) &&  (BOPAIRS2(i,1)<= 99)
       idx1 = strfind(currentstatefortest,['L0',num2str(BOPAIRS2(i,1)),'N']);
       % because states are repeated in supversiors automaton
        for z = 1:300
        TF = isempty(idx1{1,z});
        if TF == 0 
        currentstatefortest{1,z}(idx1{1,z}+4) ='T';
        end
        end
           
       else
          idx1 = strfind(currentstatefortest,['L',num2str(BOPAIRS2(i,1)),'N']);
       % because states are repeated in supversiors automaton
        for z = 1:300
        TF = isempty(idx1{1,z});
        if TF == 0 
        currentstatefortest{1,z}(idx1{1,z}+4) ='T';
        end
        end    
       end
   end
    
   if  BOPAIRS2(i,2)<10
   % label the 2nd line contingency 
           idx2 = strfind(currentstatefortest,['L','00',num2str(BOPAIRS2(i,2)),'N']);
           % because states are repeated in supversiors automaton
            for z = 1:300
            TF = isempty(idx2{1,z});
                if TF == 0 
                    currentstatefortest{1,z}(idx2{1,z}+4) ='T';
                end
            end                        
     
   else
        if (BOPAIRS2(i,2)>=10) &&  (BOPAIRS2(i,2)<= 99)
               idx2 = strfind(currentstatefortest,['L0',num2str(BOPAIRS2(i,2)),'N']);
               % because states are repeated in supversiors automaton
                for z = 1:300
                TF = isempty(idx2{1,z});
                    if TF == 0 
                        currentstatefortest{1,z}(idx2{1,z}+4) ='T';
                    end
                end    
               
        else
               idx2 = strfind(currentstatefortest,['L',num2str(BOPAIRS2(i,2)),'N']);
               % because states are repeated in supversiors automaton
                for z = 1:300
                TF = isempty(idx2{1,z});
                    if TF == 0 
                        currentstatefortest{1,z}(idx2{1,z}+4) ='T';
                    end
                end    
        end
   
   end   
   
   %% after that, we will take the tree construction process inside the 
   
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
Npointer = currentstatefortest;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
   % we  will add another for loop to change the loading conditions 
     
%% change the loading condition for the 3 loads 
    
    % assign the new load that was generated randomly  
    ps.shunt(:,2) =   [p{i,2}];  
    ps = rebalance(ps);
    pre_contingency_flows = ps.branch(:,C.br.Pf);

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% without conrtrol 
    psnocontrol = ps; 
    opt.sim.use_control = false;
    opt.sim.control_method = 'none'; 
    % outage
    br_outages = BOPAIRS2(i,:);
    % run the simulator
    fprintf('Running simulation %d',i);
    disp ('runing sumualtion without control')
    [~,relay_outages_nocontrol,MW_lost_nocontrol,p_out,busessep,flows] = dcsimsep(psnocontrol,br_outages,[],opt);
    relay_outage_nocontrol_N_2{i} = {relay_outages_nocontrol,psnocontrol.shunt(:,2),MW_lost_nocontrol}; % this line is used to record the outages for each different initial trip wasseem
%     fprintf(' Result: %.2f MW of load shedding\n',MW_lost_2(i));
%     is_blackout = dcsimsep(ps,br_outages,[],opt);
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% try again with Emergency Control control
disp('Testing DCSIMSEP with emergency control.');
opt.sim.use_control = true;
opt.sim.control_method = 'emergency_control';
pswithemergencycontrol = ps;
fprintf('Running simulation %d ',i);
opt.optimizer = 'linprog';
    
% we replaced here the cell tree object file with the tree object only, we will label the tree with required pattern inside DES Central function 
    [~,relay_outages_withemergencycontrol,MW_lost_Emergency,p_out,busessep,flows] = dcsimsep(pswithemergencycontrol,br_outages,[],opt);
    relay_outage_withemergencycontrol_N_2{i} = {relay_outages_withemergencycontrol,pswithemergencycontrol.shunt(:,2),MW_lost_Emergency};



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% try again with DES Modular control
disp('Testing DCSIMSEP with DES Modular control.');
pswithcontrol = ps;
opt.sim.control_method = 'distributed_control';
opt.sim.use_DESModular = true;
opt.optimizer = 'linprog';
opt.verbose = true;
opt.sim.nHopLoc =1;
% update supervisors from (relay_outages_withcontrol) every iteration, by updating pointers of supervisors  (pointer is nothing but the state name) 
[~,relay_outages_withcontrol,MW_lost_DES,p_out,busessep,flows] = dcsimsep(pswithcontrol,br_outages,[],opt,currentstatefortest,Npointer,pre_contingency_flows); % for DES Modular we will not need outputtreeObject variable to be passed like LLP Function
%relay_outage_withcontrol_withdata{i,j} = {relay_outages_withcontrol,pswithcontrol.shunt(:,2),[CrticalLineIndex,alcritical]};
  relay_outage_withDEScontrol__Modular_N_2{i} = {relay_outages_withcontrol,pswithcontrol.shunt(:,2),MW_lost_DES,[alIndex,alvalue,sigmaFl,uFt,pre_contingency_flows]}; % the first two variables are sorted

end


overalldataN_2 = {relay_outage_withemergencycontrol_N_2,relay_outage_nocontrol_N_2, relay_outage_withDEScontrol__Modular_N_2}; % register the data for cmparision 
% overalldataN_1 = {relay_outage_withcontrol_N_1, relay_outage_withemergencycontrol_N_1,relay_outage_nocontrol_N_1}; % register the data for cmparision 

save('overalldataN_2_30_bus_modular.mat', 'overalldataN_2', '-v7.3');
%% plot function for the N-2 cases 


y = [];
for i=1:100
    if (~isempty(overalldataN_2{1,1}{1,i})) && (~isempty(overalldataN_2{1,2}{1,i})) && (~isempty(overalldataN_2{1,3}{1,i})) 
        y = [y;overalldataN_2{1,2}{1,i}{1,3}.rebalance , (overalldataN_2{1,1}{1,i}{1,3}.rebalance+overalldataN_2{1,1}{1,i}{1,3}.control) , (overalldataN_2{1,3}{1,i}{1,3}.rebalance+overalldataN_2{1,3}{1,i}{1,3}.control)];
    end
end
bar(y);


% genearte also a box plot 

boxplot([y(:,1),y(:,2),y(:,3)])

% another way to plot all 3 cases on differnt colors with subplot
y1=12000;
subplot(1,3,1);boxplot([y(:,1)],'PlotStyle','compact','color','k'); ylim([0 y1]);subplot(1,3,2);boxplot([y(:,2)],'PlotStyle','compact','color','b'); ylim([0 y1]);subplot(1,3,3);boxplot([y(:,3)],'PlotStyle','compact','color','r'); ylim([0 y1]);
subplot(1,3,1);boxplot([y(:,1)],'PlotStyle','compact','color','k'); ylim([0 y1]);subplot(1,3,2);boxplot([y(:,2)],'PlotStyle','compact','color','r'); ylim([0 y1]);subplot(1,3,3);boxplot([y(:,3)],'PlotStyle','compact','color','b'); ylim([0 y1]);


return