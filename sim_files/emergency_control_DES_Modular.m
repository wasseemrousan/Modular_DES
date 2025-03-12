function [dPd,dPg] = emergency_control_DES_Modular(ps_agent,OptVarBus_I,opt,currentstates,ramp_limits,br_over_id) 


% initialize dPg and dPd
dPg = zeros(size(ps_agent.gen(:,1))); 
dPd = zeros(size(ps_agent.shunt(:,1)));

global pattern;
% pattern = [currentstates{1,:}];
pattern = currentstates; 
  
                %% ###############################################################################################################
                % بسم الله الرحمن الرحيم
                % 2/8/2023 : اعادة كتابة OPTIMIZATION PROBLEM IS REWRITTEN
                % TO BE MODULAR - Wasseem Al Rousan
                % ###############################################################################################################
                    
                %% 1 st : collect some data from each agent => need verification (all of the lines)
                
                C = psconstants;
                % get some information
                Np = opt.sim.Np;
                % collect some data
                nbus = size(ps_agent.bus,1);
                br_st = (ps_agent.branch(:,C.br.status) == 1);
                F = full(ps_agent.bus_i(ps_agent.branch(br_st,1)));
                T = full(ps_agent.bus_i(ps_agent.branch(br_st,2)));
                X = ps_agent.branch(br_st,C.br.X);
                G = ps_agent.bus_i(ps_agent.gen(:,1));    % generator locations
                D = ps_agent.bus_i(ps_agent.shunt(:,1));  % demand/shunt locations
                   

                % find all the generators, loads and branches in the local neighborhood
                if isempty(OptVarBus_I)
                    OptVarBus_I = 1:nbus;
                end
                
                ig = ismember(G, OptVarBus_I);
                ng = sum(ig);
                ish = ismember(D, OptVarBus_I);
                nd = sum(ish);
                genbusidx = ps_agent.gen(find(ig),1);
                shuntbusidx = ps_agent.shunt(find(ish),1);
                if ~isempty(intersect(genbusidx,shuntbusidx))
                
                end
                if ng == 0 || nd == 0
                    dPd = 0; 
                    dPg = 0;
                    return
                end
                Pg0 = ps_agent.gen(ig,C.ge.P) .* ps_agent.gen(ig,C.ge.status)./ ps_agent.baseMVA; % in per unit;
                Pd0 = ps_agent.shunt(ish,C.sh.P) .* ps_agent.shunt(ish,C.sh.factor);
                [br_id,~] = find(ismember(F,OptVarBus_I) | ismember(T,OptVarBus_I));
                nbr = length(br_id);
                f0_all = ps_agent.branch(br_st,C.br.Pf)/ps_agent.baseMVA;
                f0 = f0_all(br_id);
                fmax_all = ps_agent.branch(br_st,C.br.rateB)/ps_agent.baseMVA;
                fmax = fmax_all(br_id);
                % assign one ref bus for each island: choose a bus that is already closest to zero
                nodes = (1:nbus)'; links = [F,T];
                grid_no = find_subgraphs(nodes,links);
                n_sub = max(grid_no);
                ref_bus_i = zeros(n_sub,1);
                theta = ps_agent.bus(:,C.bu.Vang) * pi / 180;
                for grid_i = 1:n_sub
                    bus_subset = find(grid_no == grid_i);
                    theta_sub = theta(bus_subset);
                    % find the bus in this island with the smallest absolute angle
                    [~,sub_idx] = min(abs(theta_sub));
                    % Add this bus to the reference list
                    ref_bus_i(grid_i) = bus_subset(sub_idx);
                end

                % build an index into the decesion vector
                ix.dPg    = reshape(1:ng*Np, ng, Np);
                ix.dPd    = reshape(1:nd*Np, nd, Np)     + ng * Np;
                ix.dtheta = reshape(1:nbus*Np, nbus, Np) + (ng + nd) * Np;
                ix.s      = reshape(1:nbr*Np, nbr, Np)   + (ng + nd + nbus) * Np;
                ix.f      = reshape(1:nbr*Np, nbr, Np)   + (ng + nd + nbus + nbr) * Np;
                nx = (ng + nd + nbus + nbr + nbr) * Np;


                %% set up x_min, x_max
                x_min = nan(nx,1);
                x_max = nan(nx,1);
                % dPg
                x_max(ix.dPg) = 0;
                x_min(ix.dPg) = [-ps_agent.gen(ig,C.ge.ramp_rate_down); repmat(-ramp_limits(ig),Np-1,1)]; % ramping for the first time step is limited based on the agent information 
                % dPd
                x_max(ix.dPd) = 0;
                x_min(ix.dPd) = -Inf;
                % dtheta
                x_min(ix.dtheta) = -Inf; x_max(ix.dtheta) = Inf;
                x_min(ix.dtheta(ref_bus_i,:)) = 0;
                x_max(ix.dtheta(ref_bus_i,:)) = 0;
                % s
                x_min(ix.s) = 0;
                x_max(ix.s) = Inf;
                % f
                x_min(ix.f) = -Inf;
                x_max(ix.f) = Inf;
                
                % make sure the problem is not infeasible because of numerical issues
                x_min(abs(x_min) < 1e-4) = 0;
                x_max(abs(x_max) < 1e-4) = 0;
                
                %% set up the objective
                c_obj = zeros(nx,1);
                c_obj(ix.dPd) = -opt.sim.cost.load / ps_agent.baseMVA; % changed to pu in obj func, just to make it compatible with emergency_control_dec
                c_obj(ix.s) = opt.sim.cost.overload;

                 
              %% for dPd 
              Pd0 = ps_agent.shunt(ish,C.sh.P).*ps_agent.shunt(ish,C.sh.factor) ./ ps_agent.baseMVA;
              FACTOR = 1; % factor used to reduce the step down of the generators/loads
              %% Wasseem 01/19/2023 => how to find the index for label [the most loaded line in the neighborhood]
                                                  % 1st build y matrix and incidence matrix then remove
                                    % unrelated lines and nodes to the neghiborhood from it 
                                        C = psconstants; % tells me where to find my data
                                        % psorig = case300_2;  % take the MATPOWER Format of the case study 
                                        psorig = load('psorig.mat');  % take this cleaned up version of the matpower format of the IEEE 300 bus sytesm 

                                        z_ineq = fmax;
                                        
                                        As = makeIncidence(psorig.psorig.bus, psorig.psorig.branch)';
                                        % convert sparse matrix into full matrix 
                                        A = full(As)'; % A is the isdence matrix
                    %                     A(:,31) = []; 
                    %                     A(:,32) = [];
                                        nbr_all_branches = [1:size(psorig.psorig.branch,1)]';
                                        br_idcomp =  nbr_all_branches(~ismember(nbr_all_branches,br_id)); % find the branches that are not in the neighborhood
                                        index_neghoborhood_buses = full(OptVarBus_I);
                                        nbr_all_buses = [1:size(psorig.psorig.bus,1)]'; 
                                        % buil another index for generator
                                        % buses and load buses, which will
                                        % be used later after optimization
                                        % problem is solved at last step to
                                        % identify variable x
                                        % index_type_bus = zeros(1,(nd+ng));

                                        if ~isempty(intersect(genbusidx,shuntbusidx)) % Wasseem 24/1/2024
                                                index_type_bus = zeros(1,((nd+ng)-size(intersect(genbusidx,shuntbusidx),1)));
                                        else
                                                index_type_bus = zeros(1,(nd+ng));
                                        end   

                                        % we shouldn't use psorig in code sgemnt below becaus ethere some mismatch between psorig 300_ps use cases 
                                        % for h = 1:(size(index_neghoborhood_buses,1))
                                        %     if (C.psorig.bus(index_neghoborhood_buses(h),2) == 3 || psorig.psorig.bus(index_neghoborhood_buses(h),2) == 2)
                                        %         index_type_bus(h) = 1; % generator bus
                                        %     else
                                        %         index_type_bus(h) = 2; % Load bus 
                                        %     end
                                        % end

                                        
                                        for h = 1:(size(index_neghoborhood_buses,1))
                                            if ~isempty(find(ps_agent.gen(:,1)==index_neghoborhood_buses(h), 1))
                                                index_type_bus(h) = 1; % generator bus
                                            else
                                                if ~isempty(find(ps_agent.shunt(:,1)==index_neghoborhood_buses(h), 1))
                                                 index_type_bus(h) = 2; % Load bus 
                                                else
                                                 index_type_bus(h) = 4; % nither gen or Load bus 
                                                end
                                            end
                                        end


                                        index_type_bus = nonzeros(index_type_bus)';
                                        buses_not_neighborhood = nbr_all_buses(~ismember(nbr_all_buses,index_neghoborhood_buses));
                                     % do the trimming only once for Matrix
                                     % A in first iteration
                                     if size(A,1) == size(psorig.psorig.branch,1)
                                        A(br_idcomp,:)= []; 
                                        A(:,buses_not_neighborhood)=[]; % delete all un related rows and clumns with zero, the branches and buses that are not part of the neighborhood
                                        A_trimmed = A;

                                     end
                                        y =  (1./abs(psorig.psorig.branch(br_id',4))); % get the reactance of the lines, then find the admitance y for each line !
                                        % update the adimatiance of lines from vector s 
                    %                     s = psorig.psorig.branch(:,11)';
                                        s = br_st(br_id);
                    %                     s(br_idcomp)= 0;
                                        y = (y.*s);
                                        yt = diag(y); 
                                        % preliminiry variables 
                                        Atelda =  (sqrt(yt)*A_trimmed);
                                        % lets calculate H based on whole ps for now 
                                        H = makePTDF(psorig.psorig.baseMVA, psorig.psorig.bus, psorig.psorig.branch,1); % generation shift factor
                                        H_Trimmed = H; 
                                        H_Trimmed(:,buses_not_neighborhood)=[];
                                        H_Trimmed(br_idcomp,:)= [];
                                        [value,index]= max(f0); % Index for the overloaded line in th eneighbowhood wich will be identified for each agent alone
                                        % Wasseem; correction, we may want
                                        % to use the varaible br_over_id
                                        % passed from the distibusted control function

                                        % Wasseem 2/1/2024
                                        % use f0_all variable to pick the most overloaded line if br_over_id is greter than 1
                                        index = br_over_id;
                                        if size(index,1) >1
                                          [M,I] = max(abs(f0_all([index],:)));
                                          index = index(I);
                                        end
                                        f_new = H(index,:);
                                        f_new_trimmed = abs(f_new);
                                        f_new_trimmed(:,buses_not_neighborhood)=[];

                        for e = 1:size(pattern,2)
                        
                                          if index <= 9
                                                pattern{1,e} = strrep(pattern{1,e}, ['L' '00' num2str(index) 'N'], ['L' '0' num2str(index) 'T']); % Wasseem 06/17/2021 we will use pattern to search for this predicted to trip line inside each of the supervisors
                                          else
                                              if (index > 9) && (index <= 99)

                                                pattern{1,e} = strrep(pattern{1,e}, ['L0'  num2str(index) 'N'], ['L0' num2str(index) 'T']); % Wasseem 06/17/2021 we will use pattern to search for this predicted to trip line inside each of the supervisors
                                              else
                                                pattern{1,e} = strrep(pattern{1,e}, ['L'  num2str(index) 'N'], ['L' num2str(index) 'T']); % Wasseem 06/17/2021 we will use pattern to search for this predicted to trip line inside each of the supervisors
                                              
                                              end

                                          end  
                        
                        end
                        
                        
                        
                             if ps_agent.bus_id < 10 % W IS BUS NUMBER TO BE FIXED LATER ON 
                            
                            %                 disp(currentstates{1,w})
                            %                 disp('start state is')
                            %                 disp(currentstates{1,w})
                                                 supuservisorname = ['supervisor00',num2str(ps_agent.bus_id)]; 
                                                 SupInfo = load(['Bus00',num2str(ps_agent.bus_id),'supData']); 
                                                 endstates = SupervisorIterator(SupInfo,currentstates{1,ps_agent.bus_id},supuservisorname); % Wasssem 10/21/2021; use the extended supervisor instead
                            
                            
                             else
                                 if (ps_agent.bus_id > 9 &&  ps_agent.bus_id <= 99)
                                                disp('start state is')
                                                disp(currentstates{1,ps_agent.bus_id})
                                                supuservisorname = ['supervisor0',num2str(ps_agent.bus_id)]; 
                                                SupInfo = load(['Bus0',num2str(ps_agent.bus_id),'supData']);
                            
                                                endstates = SupervisorIterator(SupInfo,currentstates{1,ps_agent.bus_id},supuservisorname); % Wasssem 10/21/2021; use the extended supervisor instead
                            
                                 
                                 else
                                                disp('start state is')
                                                disp(currentstates{1,ps_agent.bus_id})
                                                supuservisorname = ['supervisor',num2str(ps_agent.bus_id)]; 
                                                SupInfo = load(['Bus',num2str(ps_agent.bus_id),'supData']);
                            
                                                endstates = SupervisorIterator(SupInfo,currentstates{1,ps_agent.bus_id},supuservisorname); % Wasssem 10/21/2021; use the extended supervisor instead
                            


                                 end
                             end
                        
                        e = cell(1,size(endstates,2));
               for r = 1:size(endstates,2)

                       e{1,r} = strfind(endstates{1,r},pattern{1,ps_agent.bus_id});
                        
                     if isempty(e{1,r})  % check over endstates if pattern (predicted line to trip is wtihin the next states of any of the supervisors) 
                               
                          disp(' NO DES Modular Control action taken')
                               
                     else

                        %% ###############################################################################################################
                        % start performing the optimization problem here % we need to make the optimation problem local and adjust the sizes of martices to suite the local neighborhood only 
                        fprintf(2,'######################################\n');
                        fprintf(2,'###DES Modular Control action taken###\n');
                        fprintf(2,'######################################\n');

                        %% Constraints
 
                                        Aineq = [(sqrt(yt)*pinv(Atelda'));f_new_trimmed]; % the decision needs to be taken must be realted to the overlaoding quantity Maybe new optimaztion problem  
                                        % a seperate 2 vectors out of x
                                        % needs to be extacted, one for
                                        % generators, the other for loads
                                        % for A_const
                                        % if ~isempty(intersect(genbusidx,shuntbusidx)) % Wasseem 24/1/2024
                                        %         A_const = ones(1,((nd+ng)-size(intersect(genbusidx,shuntbusidx),1)));
                                        % else
                                        %         A_const = ones(1,(nd+ng));
                                        % end   
                                        buscounter = size(full(OptVarBus_I),1); % lets use this counter size instead of the one below for now
   
                                        A_const = ones(1,buscounter);                                
                                           x_ig_index = 1; 
                                           x_ish_index = 1; 
                                           % if ~isempty(intersect(genbusidx,shuntbusidx)) % Wasseem 24/1/2024
                                           %   buscounter = (nd+ng)-size(intersect(genbusidx,shuntbusidx),1);
                                           % 
                                           % else
                                           %    buscounter = (nd+ng);
                                           % 
                                           % end

                                           for k =1:(buscounter)
            
                                               if (index_type_bus(k) == 1)
                                                  A_const(k)= 1;
                                                  x_ig_index = x_ig_index+1;
                                               else
                                                  A_const(k)= -1;
                                                  x_ish_index = x_ish_index+1;
                                              end
                                            end
                                        
                                        Aconst = [A_const]; 
                                        OverLoadingQuantity = abs(psorig.psorig.branch(index,8)/ps_agent.baseMVA  - abs(psorig.psorig.branch(index,C.br.Pf)/ps_agent.baseMVA)); 
                                        bconst =[0];

                                        bineq = [(z_ineq - (psorig.psorig.branch(br_id,C.br.Pf)./ps_agent.baseMVA));OverLoadingQuantity]; 
                                        % *********************************************************************
                                                     
                                        
                                    %% solve the problem
                    
                                 %% 2nd ; assume now we have subgraph that represent each agent its neighborhood 
                                 % => we need to formulate the optimization problem and get the modular supervisor info for each agent (supervisor)
                                
                                 xo = -0.0.*ones(1,(buscounter));
                                 options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',200000,'MaxIter',50000,'StepTolerance', 1e-25); 
                                            
                                  c = ones(1,(buscounter)); % c is additional cost functoin which be adjusted if needed  
                    %                   fun = @(x)(( c(1)*x(1)+c(2)*x(2) +c(3)*x(3)+c(4)*x(4)+c(5)*x(5)+c(6)*x(6)+ c(6)*x(7)+ c(8)*x(8)+ c(9)*x(9)+c(10)*x(10)+c(11)*x(11)+c(12)*x(12)+c(13)*x(13)+ c(14)*x(14)+c(15)*x(15)+c(16)*x(16)+ c(17)*x(17)+ c(18)*x(18)));
                                          
                                        syms x [1 (buscounter)]; % we will use buscounter insted of nd+ng to avoid the conflict in case we have gen and shunt at same bus

                                        for h = 1:(buscounter)                                            
                                            
                                            if index_type_bus(h) == 1
                                               c(h) = 1;
                                            end

                                        end 
                                        fun = @(x) sum(c.*x); 

                                        disp('x_min and x_max bounds are set')
                                        x_min_boundary = -ones(1,(buscounter));
                                        x_max_boundary = -0.000.*ones(1,(buscounter));
%                                         x_min_boundary(ix.dPg) = [-ps_agent.gen(ig,C.ge.ramp_rate_down); repmat(-ramp_limits(ig),Np-1,1)]; % ramping for the first time step is limited based on the agent information
                                        % set correct x_min limits for
                                        % generator buses
                                        gen_index = 1;
                                        Ld_index = 1; 
                                        gen_bus_num = zeros(1,ng);
                                        gen_bus_num_indx = 1; 
                                        for v=1:size(index_neghoborhood_buses,1)
                                            if (~isempty(find(ps_agent.gen(:,1)==index_neghoborhood_buses(v), 1)))
                                               
                                                    if gen_bus_num(gen_bus_num_indx) ==0
                                                        gen_bus_num(gen_bus_num_indx) = index_neghoborhood_buses(v);
                                                        gen_bus_num_indx = 1 + gen_bus_num_indx;  
                                                    end
                                                
                                            end
                                        end
                                        
                                        for r = 1:size(gen_bus_num,2)
                                            gen_number(r) = find(ps_agent.gen(:,1) == gen_bus_num(r));
                                        end

                                        %~isempty(find(ps_agent.gen(:,1)==index_neghoborhood_buses(h), 1))
                                        
                                        Load_bus_num = zeros(1,(nd));
                                        Ld_bus_num_indx = 1; 
                                        for v=1:size(index_neghoborhood_buses,1)
                                            if (~isempty(find(ps_agent.shunt(:,1)==index_neghoborhood_buses(v), 1))) 
                                               
                                                    if Load_bus_num(Ld_bus_num_indx) == 0
                                                        Load_bus_num(Ld_bus_num_indx) = index_neghoborhood_buses(v);
                                                        Ld_bus_num_indx = 1 + Ld_bus_num_indx;  
                                                    end
                                                
                                            end
                                        end
                                        for a=1:size(Load_bus_num,2)
                                            % new 24/1/2024 wasseem
                                            if (~isempty(find(ps_agent.shunt(:,1) == Load_bus_num(a))))
                                            Load_nmbr(a) = find(ps_agent.shunt(:,1) == Load_bus_num(a));
                                            end
                                        end
                                        
                                        sf_cap = ps_agent.capacity.sf;
                                        if( isempty(sf_cap) )
                                            sf_cap = 1;
                                        end
                                        if isempty(ps_agent.capacity.Pg)
                                            ps_agent.capacity.Pg = 0.5;
                                        end
                                        
                                        % use this function to locate if
                                        % there are any subgrids and set
                                        % load or gen bus of the other grid
                                        % to zero
                                        [sub_grids,n_sub_old] = find_subgraphs(ps_agent.bus(:,1),ps_agent.branch(br_st,1:2));
                                        % update generation capacity => line code borrowed from distributed control function
                                        ps_agent.capacity.Pg = ps_agent.gen(ig,C.ge.Pg) .* ps_agent.gen(ig,C.ge.status);
                                        pgcap_updated = ps_agent.capacity.Pg;
                                        %+++++++++++++++++++++++++++++++++++++%
                                        ge_status = ps_agent.gen(ig,C.ge.status)==1;
                                        Pg0 = ps_agent.gen(ig,C.ge.P).*ge_status./ ps_agent.baseMVA; % in per unit
                                        ramp_limits_pu = ps_agent.gen(ig,C.ge.ramp_rate_down) / ps_agent.baseMVA;
                                        Pg_max = Pg0; % < this could change
                                        Pg_min = max(0,max(Pg0-ramp_limits_pu,ps_agent.gen(ig,C.ge.Pmin).*ge_status./ps_agent.baseMVA));
                                        %+++++++++++++++++++++++++++++++++++++%
                                        % fing gen locations in the ps_agent.gen array
                                        gen_loc = find(ig==1);
                                        load_loc = find(ish==1);
                                        zax = intersect(ps_agent.gen(gen_loc,1),ps_agent.shunt(load_loc,1))';
                                        if ~isempty(zax)
                                        [efv,erf] = find(ps_agent.shunt(:,1) == zax);
                                        else 
                                          efv=0;  
                                        end
                                        load_loc = setdiff(load_loc,efv);
                                        for sz = 1:(buscounter)
                                            if index_type_bus(sz) == 1
                                                % x_min_boundary(sz) = ramp_limits((gen_number(gen_index)))*ps_agent.capacity.Pg; 
                                                % x_min_boundary(sz) = 0; 
                                                %x_min_boundary(sz)= x_max_boundary(sz)*-1; 
                                                % x_min_boundary(sz)= -ramp_limits((gen_number(gen_index)))*ps_agent.capacity.Pg*ps_agent.gen(gen_number(gen_index),2); 
                                                % x_min_boundary(sz)= -(ps_agent.capacity.Pg/ramp_limits((gen_number(gen_index)))); 
                                                % if -ps_agent.gen(gen_number(gen_index),2) >= x_min_boundary(sz) 
                                                %     x_min_boundary(sz) = -( ps_agent.gen(gen_number(gen_index),2))/ramp_limits((gen_number(gen_index)));
                                                % end
                                                
                                                x_min_boundary(sz) = -ramp_limits((gen_number(gen_index)))*ps_agent.gen(gen_index,C.ge.status)./ ps_agent.baseMVA;
                                                
                                                % ##### borrowed from the function taken from the simulator #####
                                                x_min_boundary(sz) = -ps_agent.gen(gen_loc(gen_index),C.ge.ramp_rate_down)./ ps_agent.baseMVA;
                                                x_min_boundary(sz) = min((Pg_min(find(gen_loc(gen_index)==find(ig==1))) - Pg0(find(gen_loc(gen_index)==find(ig==1))))/FACTOR,0);
                                                if ismember(ps_agent.gen(gen_loc(gen_index),1),zax) 
                                                   x_min_boundary(sz) = min((Pg_min(find(gen_loc(gen_index)==find(ig==1))) - Pg0(find(gen_loc(gen_index)==find(ig==1)))+Pd0(find(gen_loc(gen_index)==find(ig==1))))/2,0); % jus put point of intersectins to be zero for now
                                                   %x_min_boundary(sz) = 0;
                                                end

                                                if (pgcap_updated(gen_index) == 0)
                                                        x_min_boundary(sz) = 0;
                                                end
                                                if (x_min_boundary(sz) <= (-pgcap_updated(gen_index)))
                                                        x_min_boundary(sz) = (-pgcap_updated(gen_index));
                                                end
                                                if sub_grids(index_neghoborhood_buses(sz)) == 2
                                                    x_min_boundary(sz) = 0;
                                                end
                                                gen_index = gen_index + 1;
                                            else
                                                if index_type_bus(sz) == 2
                                                Load_nmbr= nonzeros(Load_nmbr)';% new 24/1/2024 wasseem
                                                if (Ld_index <= size(Load_nmbr,2))
                                                x_min_boundary(sz) = -(ps_agent.shunt((Load_nmbr(Ld_index)),C.sh.P).*sf_cap);
                                                x_min_boundary(sz)= -(Pd0(find(load_loc(Ld_index)==find(ish==1))).*sf_cap);
                                                if abs(x_min_boundary(sz)) > abs(ps_agent.shunt(find(load_loc(Ld_index)==find(ish==1)),C.sh.P).*sf_cap)
                                                  x_min_boundary(sz) = 0;
                                                end
                                                if x_min_boundary(sz) > 0
                                                  x_min_boundary(sz) = 0;
                                                end
                                                if isnan(x_min_boundary(sz))
                                                
                                                    x_min_boundary(sz) = 0;
                                                end
                                                end
                                                if ps_agent.shunt((Load_nmbr(Ld_index)),C.sh.P) < 0
                                                    x_min_boundary(sz) = 0;
                                                end
                                                if sub_grids(index_neghoborhood_buses(sz)) == 2
                                                    x_min_boundary(sz) = 0;
                                                end
                                                Ld_index = Ld_index+1;

                                                else 
                                                    x_min_boundary(sz) = 0;
                                                end
  
                                            end

                                        end
                                        
%                                         Pd_cap = ps_agent.shunt((Load_nmbr(Ld_index)),C.sh.P).*sf_cap;

                    %                     options = optimoptions("fmincon",...
                    %                         "Algorithm","interior-point",...
                    %                         "EnableFeasibilityMode",true,...
                    %                         "SubproblemAlgorithm","cg");
                    %                     x = fmincon(@(nd)fun(nd,ng),xo,Aineq,bineq,Aconst,bconst,x_min,x_max,[],options) 
           
                                         x = fmincon(fun,xo,Aineq,bineq,Aconst,bconst,x_min_boundary,x_max_boundary,[],options) 
                    
                    %                     f = ones(1,(nd+ng));
                    %                     x = linprog(f,xo,Aineq,bineq,Aconst,bconst,x_min,x_max)


                                 % initialize dPg and dPd
                    dPg = zeros(size(ps_agent.gen(:,1))); 
                    dPd = zeros(size(ps_agent.shunt(:,1)));
                               % a seperate 2 vectors out of x needs to be extacted, one for generators, the other for loads
                               x_ig = zeros(1,ng);
                               x_ish = zeros(1,nd);
                               x_ig_index = 1; 
                               x_ish_index = 1; 
                                      
                               % genbusidx = ps_agent.gen(find(ig==1),1);
                               genbusidx = gen_loc;
                               % shuntbusidx = ps_agent.shunt(find(ish==1),1);
                               shuntbusidx_new = load_loc;
                               for k =1:size(x,2)

                                   if (index_type_bus(k) == 1)
                                      x_ig(x_ig_index) = x(k);
                                      dPg(genbusidx(x_ig_index)) = x(k);
                                      x_ig_index = x_ig_index+1;
                                   end
                                   if(index_type_bus(k) == 2)
                                         x_ish(x_ish_index) = x(k);
                                         dPd(shuntbusidx_new(x_ish_index)) = x(k);
                                         x_ish_index = x_ish_index+1;
                                   
                                  end
                                
                               end
                                
                                
                                % dPd(ish==1) = x_ish;
                                % dPg(ig==1) = x_ig;
                                dPg(ig) = dPg(ig) * ps_agent.baseMVA;
                                dPd(ish) = dPd(ish) * ps_agent.baseMVA;
                                disp('assigning optimal solution to related buses in neighborhood')
%                                 if sum(dPg)==0 || sum(dPd)==0
%                                     dPg = 0;
%                                     dPd = 0;
%                                 end

%                                 break % break the for loop that goes over the endstates
                    end

              end

% end
    