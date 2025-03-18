
clc
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
opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
opt.sim.fast_ramp_mins = 1;
opt.sim.simple_redispatch = false;

%% Prepare and run the simulation for the Polish grid
% ps = case300_001_ps;

% ps = case30_ps;
ps = case118_ps;

psL = ps; % use psL for load variations purpose 

fprintf('----------------------------------------------------------\n');
disp('loading the data');
tic

ps = case118_ps;

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
MU_Loads = [ps.shunt(:,2)]'; % determine the means of the loads 

% Define the set for N-1 N-2 N-3 Outages  
BOPAIRS2 = nchoosek(1:186,2); %N-2

for j = 1:size(BOPAIRS2,1)
    % BOPAIRS2(j,:) = randperm(118,2);
    for k = 1:size(ps.shunt,1)
            psloads{1,j}(k) = normrnd(MU_Loads(k),(0.3*MU_Loads(k)))+20; 

    end
end


for i = 1:size(BOPAIRS2,1)
    
   % the current state string needs to be modifed to match the branch outage continegency 
   BOPAIRS2(i,1); % first line intentianally tripped 
   BOPAIRS2(i,2); % first line intentianally tripped

   
%% change the loading condition for the 3 loads 
    
    % assign the new load that was generated randomly  
    ps.shunt(:,2) =    [psloads{1,i}]'; 
    ps = rebalance(ps);
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% without conrtrol 
    psnocontrol = ps; 
    opt.sim.use_control = false;
    opt.sim.control_method = 'none'; 
    % outage
    br_outages = BOPAIRS2(i,:);
    % run the simulator
    fprintf('Running simulation %d  ',i);
    disp ('runing sumualtion without control')
    [~,relay_outages_nocontrol,MW_lost_nocontrol(i),p_out,busessep,flows] = dcsimsep(psnocontrol,br_outages,[],opt);
    relay_outage_nocontrol_N_2{i} = {relay_outages_nocontrol,psnocontrol.shunt(:,2)}; % this line is used to record the outages for each different initial trip wasseem
%     fprintf(' Result: %.2f MW of load shedding\n',MW_lost_2(i));
%     is_blackout = dcsimsep(ps,br_outages,[],opt);
    
end

overalldataN_2_NoControl = {relay_outage_nocontrol_N_2}; % register the data for cmparision 

% save('overalldataN_2_118_bus_newloadingcondition_NoControl.mat', 'overalldataN_2_NoControl', '-v7.3');

% check the data arrays for most vulnerable scenarios 
mostvulenrable = []; 
for i = 1:size(BOPAIRS2,1)
        if size(relay_outage_nocontrol_N_2{1,i}{1,1},1) >= 2
            mostvulenrable =[ mostvulenrable i];
        end
end
% store thew  N-2 cases 
for t= 1:size(mostvulenrable,2)
    N_2_scenario(t,:) = BOPAIRS2(mostvulenrable(t),:);
end
save ('N_2_scenario_118_bus.mat', 'N_2_scenario', '-v7.3');

%% all cases simulaTIONS BEFORE CONTROL WITH 0.7 DEVIATION 

% hold on 
% 
% for i=1:820
%     for j =1:10
%     scatter(overalldataN_2_NoControl{1,1}{i,j}{1,1}(:,1),[1:size(overalldataN_2_NoControl{1,1}{i,j}{1,1}(:,1),1)],'y'); 
%     end 
% end 
% 
% hold off


%% plot the loads for comaprasion 

% y = [];
% for i=1:10% dimension will be 163 I guess !
%     for j =1:64
%         y = [y;psloads{1,i}(1,j)];
%     end
% end
% bar(y);