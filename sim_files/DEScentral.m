
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% DEScentral function is a copy of CentralControl function, emergemcy
% control function is replaced with DEScontroller function, wehere we
% implemented our algorithm 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [ps, MW_lost, imbalance] = DEScentral(ps,sub_grids,ramp_rate,opt,pf_model,currentstates,outputtreeObject,Npointer,pre_contingency_flows,br_outages)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% these variables definitions calculate the size of the each Window N til
% reaching M 
componentsnumber = size(ps.bus,1) + size(ps.branch,1); % count the branches and the buses in the case study 
% based on the current automatas, we can detrmine the size of each window N
% for N = 1
Firstwindowsize = Npointer{1}(1);  % N =1 
% Secondwindowsize = Npointer{2}(1); % N = 2 
% Thirdwindowsize = Npointer{3}(1);
% Fourthwindowsize = Npointer{4}(1); 
% N =1 pointer 
N1pointer = Npointer{1}(1);
% Wasseem 08/11/2021, for the 118 bus system if all the buses to be taken
% in one system, we need to build the tree for one step lookahead only

% N = 2 pointer 
% N2pointer = Npointer{2}(1); 
% N = 3 pointer 
% N3pointer = Npointer{3}(1); 
% N = 4 pointer 
% N4pointer = N3pointer + Fourthwindowsize; 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


% Constants
C = psconstants;
EPS = 1e-4;
% Collect some data from the system
n = size(ps.bus,1);
m = size(ps.branch,1);
F = ps.bus_i(ps.branch(:,C.br.from));
T = ps.bus_i(ps.branch(:,C.br.to));
flow = ps.branch(:,C.br.Pf);
flow_max = ps.branch(:,opt.opf.contg_rate);
branch_st = ps.branch(:,C.br.status);
ge_status = ps.gen(:,C.ge.status);
% modified by by Wasseem 08/20/2020
% ge_status = ones(6,1);
% shall we remove this modification for the 118 bus system ? 
ge_status = ones(54,1);


Pg_max = ps.gen(:,C.ge.Pmax).*ge_status;
Pg_min = ps.gen(:,C.ge.Pmin).*ge_status;
G = ps.bus_i(ps.gen(:,1));
D = ps.bus_i(ps.shunt(:,1));
Pd0_sum = sum(ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor));
Pg0_sum = sum(ps.gen(:,C.ge.P) .* ps.gen(:,C.ge.status));
comm_status = true(n,1);

% if there are overloads in the system, try to mitigate them
if any(abs(flow)>flow_max + EPS)
    no_change = false; % things will change in ps
    % Check the mismatch (in dc only)
    if strcmp(pf_model,'dc')
        mis_old = total_P_mismatch(ps);
        if abs(mis_old)>EPS
            error('System not balanced on entry to take_control_actions');
        end
    end
    % Figure out the ramp rate
    ramp_dt = opt.sim.fast_ramp_mins * 60; 
    ramp_limits = ramp_rate*ramp_dt;
    % Find the optimal load/gen shedding
    if opt.sim.use_mpc
        [delta_Pd, delta_Pg] = solve_mpc_emergency(ps,[],ramp_limits,opt);
    else
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
global al;   
global sigmaFl;
%% 
% at some point here the markovial label file needs to be implemented. 
% we will use Markovianlabel function in a foe loop to implement backwards
% search 
% First: change the name of the case study from ps to ps1, because we wil
% do our analysis to ps1 and simulate the next two or three trips and
% inject control actions before they happen (Control actions will be
% applied to ps) 

ps1 = ps; 
% I will define the labels for N =1
% pattern = 'L01NL02NL03NL04NL05NL06TL07NL08NL09NL10NL11ND01ND02ND03NG01NG02NG03N'; 
pattern = [currentstates{1,:}]; 
% tripsN{1} = contains(outputtreeObject.Node(1: N1pointer),pattern); % for third step
% pattern = 'L01NL02NL03NL04NL05TL06TL07TL08NL09NL10NL11ND01ND02ND03NG01NG02NG03N'; 
% tripsN{2} = contains(outputtreeObject.Node(N1pointer: (N2pointer-1)),pattern);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

[CptMatrix,al,sigmaFl] = Markovianlabelforcontrol(outputtreeObject,ps,pre_contingency_flows,br_outages); % this function is a copy of Markovianlabel

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% this part here to be used later to label illegal states for the next look ahead windows
% Wasseem 08/11/2021, for the 118 bus system if all the buses to be taken in one system, we need to build the tree for one step lookahead only
     
[alcriticalforlabel,Indexforlabel] = Markovianlabel(outputtreeObject,ps,CptMatrix,al); % use the raltions described in Markovian transition model paper  
%% trip the crtical line 
% disp('line is')
% disp(Indexforlabel)

ps1.branch(Indexforlabel,C.br.status) = 0;  
% get l critical info and 

%% Label the output tree file, label the illegal state for next step.
% use strrep; Find and replace substrings
if Indexforlabel < 10
    pattern = strrep(pattern, ['L' '00' num2str(Indexforlabel) 'N'], ['L' '00' num2str(Indexforlabel) 'T']);
else
    if (10 <= Indexforlabel) && (Indexforlabel< 100)

        pattern = strrep(pattern, ['L' '0' num2str(Indexforlabel) 'N'], ['L' '0' num2str(Indexforlabel) 'T']);
    
    else
        pattern = strrep(pattern, ['L' num2str(Indexforlabel) 'N'], ['L' num2str(Indexforlabel) 'T']);

    end
end
% wasseem 08/11/2021 , the h == 2 is removed from this part here beacuse we
% consider one step lookaheaindow only 
tripsN{1} = contains(outputtreeObject.Node(1: N1pointer),pattern); % for third step
   
% al to be modified for later steps to include the more than one step
% (sequnbtial tripping) 
% this may need to perform a for loop inside this function for 3 - 4 steps
% for example , and we will get a string for several sequential critical
% trips 


labeledtree = {outputtreeObject, tripsN};  
% save labeledtree; 

%%
% global CrticalLineIndex;
% global alcritical;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%


[delta_Pd,delta_Pg] = takecontrolactions(labeledtree,currentstates,ps,Npointer,al,sigmaFl); % this to be implemented as acallback function


%% we will need to import some functions from Emergency control or Central Control here to implement the Control Actions  
        
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    end
    % If emergency control says that we should do something:
    if any(abs(delta_Pd)>EPS)
        if opt.verbose, print_load_shedding_results(ps,delta_Pg,delta_Pd); end
        % Compute the new amount of generation
        Pg = ps.gen(:,C.ge.P) + delta_Pg;
        disp(' Pg is :')
        disp(Pg);
        % check if Pg is in its bounds if needed
        if opt.pf.check_Pg
            if any( Pg>Pg_max+EPS | Pg<Pg_min-EPS )
%                 disp ('Pg_max+EPS is ')
%                 disp(Pg_max+EPS)
%                 disp('Pg_min-EPS is')
%                 disp(Pg_min-EPS)
                error('Pg is out of bounds'); 
            end
        end
        % Compute the new load factor
        delta_lf = delta_Pd./ps.shunt(:,C.sh.P); % in IEEE 30 BUS case we have shunt nodes with zero P load, this will lead to make delta_lf = -inf
		%++++++++++++++%
% 		for i = 1:size(delta_lf,1)
% 			if delta_lf(i) == -inf;
% 				delta_lf(i) = -0.000001;
% 			end
% 		end
		%++++++++++++++%
		
        disp('delta_lf');
        disp(delta_lf);
        delta_lf(isnan(delta_lf)) = 0;
        % make sure ps.shunt(:,C.sh.factor) are all ones; Wasseem
        % 06/18/2020, 
        ps.shunt(:,C.sh.factor) = ones(size(ps.shunt(:,C.sh.factor),1),1); 
%         ps.shunt(:,C.sh.factor) = ones(24,1);
%         lf_new = ps.shunt(:,C.sh.factor) + delta_lf; % there is a problem in this part : ps.shunt(:,C.sh.factor) return sometimes [ 0 1 1]', in the priginal system it is  [1 1 1]'
        lf_new = ones(1,64)' + delta_lf; % added by wasseem 05/10/2020 , modfied at 08/15/2021 to be 64 for the 118 BUS system
%         this part was added for the ieee6ww bus system 
%         disp('ps.shunt(:,C.sh.factor)')
%         disp(ps.shunt(:,C.sh.factor))
%         disp('lf_new with DES control')
%         disp(lf_new)
        if any(lf_new < -EPS | lf_new > 1)
            error('Load factors are out of bounds.')
        end
        % Implement the results
        ps.gen(:,C.ge.P) = Pg; % implement Pg
        ps.shunt(:,C.sh.factor) = lf_new; % implement shunt factor % wasseem 
        % record new Pg and Pd
        Pd_sum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
        Pg_sum = sum(ps.gen(:,C.ge.P).*ps.gen(:,C.ge.status));
        % Get the imbalance caused by control 
        imbalance = abs( (Pg_sum-Pg0_sum) - (Pd_sum-Pd0_sum) );
        if imbalance < EPS
            imbalance = 0;
        else
            % the system needs a rebalance
            ps = rebalance(ps,sub_grids,ramp_limits,opt,[],pf_model);
            Pd_sum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
        end
        MW_lost.control = Pd0_sum - Pd_sum;
        % run power flow
        if strcmp(pf_model,'dc')
            ps = dcpf(ps,sub_grids);
            MW_lost_total = MW_lost.control;
        elseif strcmp(pf_model,'ac')
            Pd_sum_before = Pd_sum;
            opt.pf.CascadeControl = true;
            ps = acpf(ps,[],opt);
            Pd_sum = sum(ps.shunt(:,C.sh.P).*ps.shunt(:,C.sh.factor));
            MW_lost.voltage_collapse = Pd_sum_before - Pd_sum;
            MW_lost_total = MW_lost.control + MW_lost.voltage_collapse;
        else
            error('Unknown power flow model.')
        end
        % check if Pg is in its bounds if needed
        if opt.pf.check_Pg
            Pg = ps.gen(:,C.ge.P);
            if any( Pg>Pg_max+EPS | Pg<Pg_min-EPS )
                error('Pg is out of bounds'); 
            end
        end
        % error check
        if abs(MW_lost_total - (Pd0_sum - Pd_sum)) > EPS
            error('Something is wrong!')
        end
    else
        no_change = true;
    end
else
    no_change = true;
end

if no_change
    MW_lost.control = 0;
    imbalance = 0;
    if strcmp(pf_model,'ac')
        MW_lost.voltage_collapse = 0;
    end
end




