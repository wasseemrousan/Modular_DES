%% this part is for calculating the control actions using LLP policy 
%tak
% econtrolactions function will be a replica of the emergeny_control function
%dcsemsep library, am also working on a method to use tree search and leabeling method   

% after buidling the tree object file, the next step is to label the tree
% nodes that reflects the components that to be tripped, in order to do
% this we need to do 2 steps: 

% 1st : search the tree for the a specific string or patttern and obtain
% the index of this node , this is already done in buildtreeobjectfile part
% , we have now 4 variables : indexN1 indexN2 indexN3 and indexN4

% process these illegal states and take control actions  

% for object file t
function [delta_Pd,delta_Pg] = takecontrolactions(labeledtree,currentstates,ps,Npointer,al,sigmaFl) % this to be implemented as acallback function

% al to be modified for later steps to include the more than one step
% (sequential tripping) 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% these variables definitions calculate the size of the each Window N til
% reaching M 
% based on the current automatas, we can detrmine the size of each window N
% for N = 1
Firstwindowsize = Npointer{1}(1);  % N =1 
% Secondwindowsize = Npointer{2}(1); % N = 2 
% Thirdwindowsize = Npointer{3}(1);
% Fourthwindowsize = componentsnumber^4; 

% N =1 pointer 
N1pointer = Npointer{1}(1); 
% N = 2 pointer 
% N2pointer = Npointer{2}(1); 
% N = 3 pointer 
% N3pointer = Npointer{3}(1); 
% N = 4 pointer 
% N4pointer = N3pointer  + Fourthwindowsize;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% to cacualte the load shedding, a method will be used based on PTDF
% calculations. 
% The makePTDF function returns the DC PTDF matrix for a given choice of slack.
% The matrix is nbr × nb, where nbr is the number of branches and nb is the number
% of buses.

% H = makePTDF(ps.baseMVA, ps.bus, ps.branch,1); 
% % get the crticial line from the labeled tree object file 
% % labeledtree 
% % the rows presents the lines, so the critical line will be the row number
% % choose the most sensitive load in that row ? 
% 
% % employ the H matrix calculation in the load shedding or generation
% % redispatch
% 
% f =  H(Index,ps.shunt(:,1)'); % we will take the load busses  
global alIndex;
global alvalue;
global Aineq;
global yt;
global Atelda; 

[alvalue,alIndex]= sort(al); % sort the whole al vector from lower to higher
C = psconstants;




psorig = case118; % take the MATPOWER Format of the case study 
psorig.branch(:,11) = ps.branch(:,11); % update branch status for the case
As = makeIncidence(psorig.bus, psorig.branch)';
% convert sparse matrix into full matrix 
A = full(As)'; % A is the isdence matrix
y =  (1./ps.branch(:,4))'; % get the reactance of the lines, then find the admitance y for each line !
% update the adimatiance of lines from vector s 
s = ps.branch(:,11)'; 
y = (y.*s);
yt = diag(y); 
% preliminiry variables 
Atelda =  (sqrt(yt)*A);
H = makePTDF(psorig.baseMVA, psorig.bus, psorig.branch,69); % generation shift factor 
% H = makePTDF(ps.baseMVA, ps.bus, ps.branch,1); % generation shift factor 

% we will use the alindex variable and take the first two most critical
% lines and take the average of them, Wasseem 06/30/2020
% f = zeros(2,30);
% f(1,:) = H(alIndex(1)); 
% f(2,:) = H(alIndex(2));
% favg = mean(f);
% this vector gets the loading (generation or loads) rating for all the busses
psforconstraints = case118_ps; 
% LoadRatings = [psforconstraints.gen(1,2),psforconstraints.gen(2,2),ps.shunt(1,2),ps.shunt(2,2),ps.shunt(3,2),ps.shunt(4,2),ps.shunt(5,2),ps.shunt(6,2),ps.shunt(7,2),ps.shunt(8,2),ps.shunt(9,2),ps.shunt(10,2),psforconstraints.gen(6,2),ps.shunt(11,2),ps.shunt(12,2),ps.shunt(13,2),ps.shunt(14,2),ps.shunt(15,2),ps.shunt(16,2),ps.shunt(17,2),ps.shunt(18,2),psforconstraints.gen(3,2),psforconstraints.gen(5,2),ps.shunt(19,2),ps.shunt(20,2),ps.shunt(21,2),psforconstraints.gen(4,2),ps.shunt(22,2),ps.shunt(23,2),ps.shunt(24,2)]';
LoadRatings = psorig.bus(:,3).*0;
%% Wasseem 08/15/2021 set correct load ratings for the load buses 

for g = 1: size(ps.shunt(:,1),1)
    
     
        LoadRatings (psorig.bus(ps.shunt(g,1),1)) = ps.shunt(g,2);
        
     
end 
% f(1,:) = H(alIndex(1),[ps.shunt(1,1),ps.shunt(2,1),ps.shunt(3,1),ps.shunt(4,1),ps.shunt(5,1),ps.shunt(6,1),ps.shunt(7,1),ps.shunt(8,1),ps.shunt(9,1),ps.shunt(10,1),ps.shunt(11,1),ps.shunt(12,1),ps.shunt(13,1),ps.shunt(14,1),ps.shunt(15,1),ps.shunt(16,1),ps.shunt(17,1),ps.shunt(18,1),ps.shunt(19,1),ps.shunt(20,1),ps.shunt(21,1),ps.shunt(22,1),ps.shunt(23,1),ps.shunt(24,1)]);
% f(2,:) = H(alIndex(2),[ps.shunt(1,1),ps.shunt(2,1),ps.shunt(3,1),ps.shunt(4,1),ps.shunt(5,1),ps.shunt(6,1),ps.shunt(7,1),ps.shunt(8,1),ps.shunt(9,1),ps.shunt(10,1),ps.shunt(11,1),ps.shunt(12,1),ps.shunt(13,1),ps.shunt(14,1),ps.shunt(15,1),ps.shunt(16,1),ps.shunt(17,1),ps.shunt(18,1),ps.shunt(19,1),ps.shunt(20,1),ps.shunt(21,1),ps.shunt(22,1),ps.shunt(23,1),ps.shunt(24,1)]);
% 
% favg = mean(f); % we may change the from the average to the summation 
% old statement: 
% LoadsIndexVector = [ps.shunt(1,1),ps.shunt(2,1),ps.shunt(3,1),ps.shunt(4,1),ps.shunt(5,1),ps.shunt(6,1),ps.shunt(7,1),ps.shunt(8,1),ps.shunt(9,1),ps.shunt(10,1),ps.shunt(11,1),ps.shunt(12,1),ps.gen(13,1),ps.shunt(14,1),ps.shunt(15,1),ps.shunt(16,1),ps.shunt(17,1),ps.shunt(18,1),ps.shunt(19,1),ps.shunt(20,1),ps.shunt(21,1),ps.shunt(22,1),ps.shunt(23,1),ps.shunt(24,1)];
% disp('PTDF factors are :')
% H(alIndex(1),:)

f = ones(1,118);
[value,BusIndex] = sort(H(alIndex(1),:)); 



[Maxvalue, Maxindex] = max(abs(H(alIndex(1),:)));


OverLoadingQuantity = (alvalue(1)*sigmaFl(alIndex(1))); 
% disp('overloading quanuity from al and sugma is:')
% OverLoadingQuantity

OverLoadingQuantity = ps.branch(alIndex(1),8)/100  - abs(ps.branch(alIndex(1),C.br.Pf)/100); % Wasseem 08/18/2020 , another way to look into the problrm !
% disp('New OverLoadingQuantity')
% disp(OverLoadingQuantity)
f = H(alIndex(1),:);%f is redifined to be PTDF factors


ub = ((0.999)).*(LoadRatings./100); % realte the upper bound to the amount of over loading may help !

lb = 0.0000000001.*(LoadRatings./100); % 07/18/2020, new formualtion of load rating


% make upper and lower bounds for the egenrators types zero then change
% them later 
 for y = 1:size(ps.gen(:,1),1)
     
      if  ps.gen(y,2) == 0
          ub(ps.bus(ps.gen(y,1))) = 0; 
          lb(ps.bus(ps.gen(y,1))) = 0; 

      end 
     
     
     
 end
% ub(2) = 0.2*(ps.shunt(1,2)/100); 
% ub(23) = 0.2*(ps.shunt(20,2)/100); 
% define some variables for optimization
ge_status = ps.gen(:,C.ge.status);
Pg_max = ps.gen(:,C.ge.Pmax).*ge_status; % copied from DES central function 
Pg_min = ps.gen(:,C.ge.Pmin).*ge_status; % copied from DES central function
% % upper bounds 
ub(10) = ((0.0))*abs(ps.gen(5,2)/100 - ps.gen(1,9)/100); 
ub(12) = ((0.0))*abs(ps.gen(6,2)/100 - ps.gen(1,9)/100); 
ub(25) = ((0.0))*abs(ps.gen(11,2)/100 - ps.gen(1,9)/100); 
ub(26) = ((0.0))*abs(ps.gen(12,2)/100 - ps.gen(1,9)/100); 
ub(31) = ((0.0))*abs(ps.gen(14,2)/100 - ps.gen(1,9)/100); 
ub(46) = ((0.0))*abs(ps.gen(20,2)/100 - ps.gen(1,9)/100); 
ub(49) = ((0.0))*abs(ps.gen(21,2)/100 - ps.gen(1,9)/100); 
ub(54) = ((0.0))*abs(ps.gen(22,2)/100 - ps.gen(1,9)/100);
ub(59) = ((0.0))*abs(ps.gen(25,2)/100 - ps.gen(1,9)/100); 
ub(61) = ((0.0))*abs(ps.gen(26,2)/100 - ps.gen(1,9)/100); 
ub(65) = ((0.0))*abs(ps.gen(28,2)/100 - ps.gen(1,9)/100);
ub(66) = ((0.0))*abs(ps.gen(29,2)/100 - ps.gen(1,9)/100);
ub(69) = ((0.0))*abs(ps.gen(30,2)/100 - ps.gen(1,9)/100);
ub(80) = ((0.0))*abs(ps.gen(37,2)/100 - ps.gen(1,9)/100); 
ub(87) = ((0.0))*abs(ps.gen(39,2)/100 - ps.gen(1,9)/100); 
ub(89) = ((0.0))*abs(ps.gen(40,2)/100 - ps.gen(1,9)/100); 
ub(100) = ((0.0))*abs(ps.gen(45,2)/100 - ps.gen(1,9)/100); 
ub(103) = ((0.0))*abs(ps.gen(46,2)/100 - ps.gen(1,9)/100); 
ub(111) = ((0.0))*abs(ps.gen(51,2)/100 - ps.gen(1,9)/100); 
% % lower bounds 
lb(10) = -((0.99))*abs(ps.gen(5,2)/100 - 0); 
lb(12) = -((0.99))*abs(ps.gen(6,2)/100 - 0);  
lb(25) = -((0.99))*abs(ps.gen(11,2)/100 - 0);  
lb(26) = -((0.99))*abs(ps.gen(12,2)/100 - 0);  
lb(31) = -((0.99))*abs(ps.gen(14,2)/100 - 0);  
lb(46) = -((0.99))*abs(ps.gen(20,2)/100 - 0);  
lb(49) = -((0.99))*abs(ps.gen(21,2)/100 - 0);  
lb(54) = -((0.99))*abs(ps.gen(22,2)/100 - 0); 
lb(59) = -((0.99))*abs(ps.gen(25,2)/100 - 0); 
lb(61) = -((0.99))*abs(ps.gen(26,2)/100 - 0); 
lb(65) = -((0.99))*abs(ps.gen(28,2)/100 - 0); 
lb(66) = -((0.99))*abs(ps.gen(29,2)/100 - 0); 
lb(69) = -((0.99))*abs(ps.gen(30,2)/100 - 0); 
lb(80) = -((0.99))*abs(ps.gen(37,2)/100 - 0);
lb(87) = -((0.99))*abs(ps.gen(39,2)/100 - 0); 
lb(89) = -((0.99))*abs(ps.gen(40,2)/100 - 0); 
lb(100) = -((0.99))*abs(ps.gen(45,2)/100 - 0); 
lb(103) = -((0.99))*abs(ps.gen(46,2)/100 - 0); 
lb(111) = -((0.99))*abs(ps.gen(51,2)/100 - 0); 

% define function to reduce loads shedding over all loads buses
fun = @(x)(( x(2)+x(3) + x(5)+x(7)+ x(9) + x(11)+ x(13)+ x(14)+ x(16) + x(17) +x(20)+ x(21)+ x(22) + x(23)+ x(28)+x(29)+ x(30)+ x(33)+x(35)+ x(37) + x(38)+ x(39) + x(41)+ x(43)...
+ x(44)+x(45) + x(47)+x(48)+ x(50) + x(51)+ x(52)+ x(53) + x(57) +x(58)+ x(60)+ x(63) + x(64)+ x(67)+x(68)+ x(71)+ x(75)+x(78)+ x(79) + x(81)+ x(82) + x(83)+ x(84)...
+ x(86)+x(88) + x(93)+x(94)+ x(95) + x(96)+ x(97)+ x(98)+ x(101) + x(102) +x(106)+ x(108)+ x(109) + x(114)+ x(115)+x(117)+ x(118) ));
% fun = @(x)(( x(1)+x(2)+x(3)+x(4) + x(5)+x(6)+ x(7) + x(8)+ x(9)+ x(10)+ x(11) + x(12) +x(13) +x(14)+ x(15)+ x(16) + x(17)+ x(18)+x(19)+ x(20)+ x(21)+x(22) +x(23)+x(24)+ x(25) + x(26)+x(27)+ x(28) + x(29)+ x(30)));

% fun = @(x)( f(1)*x(1)+f(2)*x(2)+f(3)*x(3)+f(4)*x(4) + f(5)*x(5)+f(6)*x(6)+f(7)*x(7)+ f(8)*x(8)+ f(9)*x(9)+ f(10)*x(10)+ f(11)*x(11)+ f(12)*x(12)+f(13)*x(13) + f(14)*x(14)+ f(15)*x(15)+ f(16)*x(16) +f(17)*x(17)+f(18)*x(18)+f(19)*x(19)+ f(20)*x(20)+ f(21)*x(21)+f(22)*x(22)+f(23)*x(23)+f(24)*x(24)+f(25)*x(25)+f(26)*x(26)+f(27)*x(27)+f(28)*x(28)+ f(29)*x(29)+f(30)*x(30));

xo = 0.001.*ones(1,118);
% xo = 0.15*(LoadRatings./100);
% [m,i] = max(abs(f)); % old statement 
% [m,i] = max(favg);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% constraints part: 07/18/2020 

% Aconst = [H(alIndex(1),:);H(alIndex(2),:);H(alIndex(3),:)]; % constraints A Matrix
% bconst = [1.2*(alvalue(1)*sigmaFl(alIndex(1)));1.2*(alvalue(2)*sigmaFl(alIndex(2)));1.2*(alvalue(3)*sigmaFl(alIndex(3)))]; % constraints b Matrix

% Aconst = [H(alIndex(1),:);H(alIndex(2),:);H(alIndex(3),:)]; % constraints A Matrix
% bconst = [(alvalue(1)*sigmaFl(alIndex(1)));(alvalue(2)*sigmaFl(alIndex(2)));(alvalue(3)*sigmaFl(alIndex(3)))];

Aconst = ones(1,118); % ones 
bconst = [0]; 
z_ineq = 0.9999.*(ps.branch(:,8))./100; 

% constraints will be related to line flows 
% realte xi index of buses with the line flows  
% do some flow calculations : s = ps.branch(:,11)'; 
% A: incidence matrix 

% P = ; % power injections at each bus
% PlFlow = sqrt(yt)*pinv(Atelda');
% DeltaP = PlFlow - ps.branch(:,8);

% [v,i] = sort((H(alIndex(2),:))); % if there is an error in the prediction algorithm it usally the second one that is correct 
% f(i(1)) = (0.2*v(1));
    
% for t =1:30
%     if (ps.bus(t,2)== 1)
%         
%         f(t) = abs(H(alIndex(1),t));
%         
%     else 
%         f(t) = -abs(H(alIndex(1),t));
%         
%     end  
% end 
    if (((value(1)))- (value(2))) <= -0.99 
%         if (ps.bus(BusIndex(1),2) == 2)||(ps.bus(BusIndex(1),2) == 3)
% %           index = ps.gen(1,BusIndex(1)); 
% %            
%             ind = find(ps.gen == BusIndex(1)); 
%             ub(BusIndex(1)) = ((0.999))*abs(ps.gen(ind,2)/100 - ps.gen(ind,9)/100);
% 
%             lb(BusIndex(1))=0;
%             disp('new ub is applied')
             f = ones(1,118);        
             f(BusIndex(1)) = 0; 
             disp(' this is where PTDF values are computed in wrong way')            
%         end
        
    end
       if (((value(118)))- (value(117))) >= 0.99
    %         if (ps.bus(BusIndex(1),2) == 2)||(ps.bus(BusIndex(1),2) == 3)
    % %           index = ps.gen(1,BusIndex(1)); 
    % %            
    %             ind = find(ps.gen == BusIndex(1)); 
    %             ub(BusIndex(1)) = ((0.999))*abs(ps.gen(ind,2)/100 - ps.gen(ind,9)/100);
    % 
    %             lb(BusIndex(1))=0;
    %             disp('new ub is applied')
                 f = -1.*ones(1,118);        
                 f(BusIndex(1)) = 0; 
                 disp(' this is where PTDF values are computed in wrong way')    
%                 ub(1) = ((0.5))*abs(ps.gen(1,2)/100 - ps.gen(1,9)/100); 
%                 ub(2) = ((0.5))*abs(ps.gen(2,2)/100 - ps.gen(2,9)/100); 
%                 ub(22) = ((0.5))*abs(ps.gen(3,2)/100 - ps.gen(3,9)/100); 
%                 ub(27) = ((0.5))*abs(ps.gen(4,2)/100 - ps.gen(4,9)/100); 
%                 ub(23) = ((0.5))*abs(ps.gen(5,2)/100 - ps.gen(5,9)/100); 
%                 ub(13) = ((0.5))*abs(ps.gen(6,2)/100 - ps.gen(6,9)/100);  
    %         end

        end
   
Aineq = [(sqrt(yt)*pinv(Atelda',1e-4));(f)];


% for n = 1:size(Aineq,1)
%     for t = 1:size(Aineq,2)
%         if (Aineq(n,t) > -0.0001)  && (Aineq(n,t) < 0)
%             Aineq(n,t) = 0;
%         end
%         
%         if (Aineq(n,t) < 0.0001) && (Aineq(n,t) > 0)
%             Aineq(n,t) = 0;
%         end
%         
%     end
% end
for n = 1:118
        if (Aineq(187,n) > -0.0001)  && (Aineq(187,n) <= 0)
            Aineq(187,n) = 1;
        end
        
        if (Aineq(187,n) < 0.0001) && (Aineq(187,n) >= 0)
            Aineq(187,n) = 1;
        end
        
end

% disp('check detrminant of (sqrt(yt)*pinv Atelda ')
% disp(det((sqrt(yt)*pinv(Atelda',1e-4))))
% 
% disp('check symbolic of detrminant of (sqrt(yt)*pinv Atelda ')
% disp(sym(det((sqrt(yt)*pinv(Atelda',1e-4)))))

% Aineq = [(f);ones(1,30)];

% bineq = [[sign(ps.branch(alIndex(1),C.br.Pf))*(alvalue(1)*sigmaFl(alIndex(1)));zeros(40,1)]+(ps.branch(:,C.br.Pf))./100;(alvalue(1)*sigmaFl(alIndex(1)))]; 
% z_ineq = zeros(41,1);
% z_ineq = al.*sigmaFl;
% z_ineq(alIndex(1)) = abs((alvalue(1)*sigmaFl(alIndex(1)))); 
% z_ineq(alIndex(2)) = abs(alvalue(2)*sigmaFl(alIndex(2))); 
% if abs(alvalue(1)*sigmaFl(alIndex(1))) >= 0.01
%     OverLoadingQuantity = abs(alvalue(1)*sigmaFl(alIndex(1))); 
% else
%     OverLoadingQuantity = 0.01; 
% end
OverLoadingQuantity = ps.branch(alIndex(1),8)/100  - abs(ps.branch(alIndex(1),C.br.Pf)/100); 

bineq = [(z_ineq - abs(ps.branch(:,C.br.Pf)./100));(OverLoadingQuantity)]; % I tried to add this to account for prediction error 2*abs(alvalue(1)*sigmaFl(alIndex(1)))+0.1*abs(alvalue(2)*sigmaFl(alIndex(2)))

    
% bineq = [(OverLoadingQuantity);10*OverLoadingQuantity]; % I tried to add this to account for prediction error 2*abs(alvalue(1)*sigmaFl(alIndex(1)))+0.1*abs(alvalue(2)*sigmaFl(alIndex(2)))




 

%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    % Wasseem 08/22/2021 : Do we need to make offline simulation to check
    % if the predicted trip line in next step is going to cuase further trips ? 
    % simulate without conrtrol 
    
    % set some options
        opt = psoptions;
        opt.verbose = false; % set this to false if you don't want stuff on the command line
        % Stopping criterion: (set to zero to simulate a complete cascade)
        opt.sim.stop_threshold = 0.00; % the fraction of nodes, at which to declare a major separation
        opt.sim.fast_ramp_mins = 1;
        opt.sim.simple_redispatch = false;
        psnocontrol = ps; 
        opt.sim.use_control = false;
        opt.sim.control_method = 'none'; 
        % Apply outage, this outage here will be the predicted to trip line
        br_outages = alIndex(1); 
        % run the simulator
        
        br_st = ps.branch(:,C.br.status)~=0;
        % check to make sure that the base case is load balanced

        [sub_grids,n_sub_old] = find_subgraphs(ps.bus(:,1),ps.branch(br_st,1:2));
%         save sub_grids; 
        save('sub_grids', 'sub_grids');
        
%         relay_outages_nocontrol_check = zeros(0);
        
        if   n_sub_old>1   % don't do optimization in this case , else do optimization as typical
        disp ('the power system has more than one island ; # of islands are: ') 
        disp (n_sub_old)
        
        % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  %
        
        %% edit this code so that it operates separately on each sub-grid
                        
        % @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ %
            x = zeros (1,118);
        % Wasseem 10/03/2021 ; this part was before using the MonteCarlo
        % Sim; for some cases it throughs an error ; Error using dcsimsep (line 94) The base case has more than one island
        delta_Pd_pu = [x(2),x(3) , x(5),x(7), x(9) , x(11), x(13), x(14), x(16) , x(17) ,x(20), x(21), x(22) , x(23), x(28),x(29), x(30), x(33),x(35), x(37) , x(38), x(39) , x(41), x(43)...
                                            , x(44),x(45) , x(47),x(48), x(50) , x(51), x(52), x(53) , x(57) ,x(58), x(60), x(63) , x(64), x(67),x(68), x(71), x(75),x(78), x(79) , x(81), x(82) , x(83), x(84)...
                                            , x(86),x(88) , x(93),x(94), x(95) , x(96), x(97), x(98), x(101) , x(102) ,x(106), x(108), x(109) , x(114), x(115),x(117), x(118)]'; % apply the optimal solution that was found 
                    %        Pd = Pd0 + delta_Pd_pu;
                             delta_Pd = delta_Pd_pu.*100; % Wasseem 07/18/2020, the base MVA is 100 , 
                             delta_Pg_pu = -1.*[x(1), x(4), x(6),	x(8),	x(10),	x(12),	x(15),	x(18),	x(19),	x(24), x(25),	x(26),	x(27),	x(31),	x(32),	x(34),	x(36),	x(40),	x(42),	x(46),	x(49) ...
                                            , x(54),	x(55),	x(56),	x(59),	x(61),	x(62),	x(65),	x(66),	x(69),	x(70),	x(72),	x(73),	x(74),	x(76),	x(77),	x(80),	x(85),	x(87),	x(89),	x(90),	x(91),	x(92) ...	
                                            , x(99),	x(100),	x(103),	x(104),	x(105),	x(107),	x(110),	x(111),	x(112),	x(113),	x(116)]'; 
                             delta_Pg = delta_Pg_pu.*100;
        
        else
        
%             psnocontrol.branch(br_outages,C.br.status) = 0;
%             br_st_nocontrol = psnocontrol.branch(:,C.br.status)~=0;
%             [sub_gridsold,n_sub_old_nocontrol] = find_subgraphs(psnocontrol.bus(:,1),psnocontrol.branch(br_st_nocontrol,1:2));
            
%             if n_sub_old_nocontrol == 1
                disp ('runing sumualtion without control')
                [~,relay_outages_nocontrol_check,MW_lost_nocontrol_check,p_out_check,busessep_check,flows_check] = dcsimsep(psnocontrol,br_outages,[],opt);
                if isempty(relay_outages_nocontrol_check)
                disp ('this continegncy will not escalate to failure cascade')    
                    
                x = zeros (1,118);
                delta_Pd_pu = [x(2),x(3) , x(5),x(7), x(9) , x(11), x(13), x(14), x(16) , x(17) ,x(20), x(21), x(22) , x(23), x(28),x(29), x(30), x(33),x(35), x(37) , x(38), x(39) , x(41), x(43)...
                               , x(44),x(45) , x(47),x(48), x(50) , x(51), x(52), x(53) , x(57) ,x(58), x(60), x(63) , x(64), x(67),x(68), x(71), x(75),x(78), x(79) , x(81), x(82) , x(83), x(84)...
                               , x(86),x(88) , x(93),x(94), x(95) , x(96), x(97), x(98), x(101) , x(102) ,x(106), x(108), x(109) , x(114), x(115),x(117), x(118)]'; % apply the optimal solution that was found 
                    %        Pd = Pd0 + delta_Pd_pu;
                delta_Pd = delta_Pd_pu.*100; % Wasseem 07/18/2020, the base MVA is 100 , 
                delta_Pg_pu = -1.*[x(1), x(4), x(6),	x(8),	x(10),	x(12),	x(15),	x(18),	x(19),	x(24), x(25),	x(26),	x(27),	x(31),	x(32),	x(34),	x(36),	x(40),	x(42),	x(46),	x(49) ...
                                   , x(54),	x(55),	x(56),	x(59),	x(61),	x(62),	x(65),	x(66),	x(69),	x(70),	x(72),	x(73),	x(74),	x(76),	x(77),	x(80),	x(85),	x(87),	x(89),	x(90),	x(91),	x(92) ...	
                                   , x(99),	x(100),	x(103),	x(104),	x(105),	x(107),	x(110),	x(111),	x(112),	x(113),	x(116)]'; 
                delta_Pg = delta_Pg_pu.*100;
                else
                    
                    
                    
                    
                   %% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++% 
       
    
                    % fmincon_option=optimset( 'Algorithm', 'interior-point','MaxFunEvals', 200000,'MaxIter',50000,'TolFun',1e-6,'TolX',1e-12);
                        options = optimoptions(@fmincon,'Algorithm', 'sqp','MaxFunctionEvaluations',200000,'MaxIter',50000,'StepTolerance', 1e-25); 
                        
                        x = fmincon(fun,xo,Aineq,bineq,Aconst,bconst,lb,ub,[],options)
                        % ##############################################################%

                         x = -x; 


                    %%
                    % labeledtree has two components, the 1st one is the tree object file, and
                    % the 2nd one is the nodes that contains the illegal states(Pattern) at each N  
                    % for N =1 
                    for i = 1:Firstwindowsize % i =1 is the root of the tree





                        if labeledtree{1,2}{1,1}(i) == 1 
                    %      disp('hiN1') %Wasseem 04/26/2020 I commnetd this so the screen will
                    

                             delta_Pd_pu = [x(2),x(3) , x(5),x(7), x(9) , x(11), x(13), x(14), x(16) , x(17) ,x(20), x(21), x(22) , x(23), x(28),x(29), x(30), x(33),x(35), x(37) , x(38), x(39) , x(41), x(43)...
                                            , x(44),x(45) , x(47),x(48), x(50) , x(51), x(52), x(53) , x(57) ,x(58), x(60), x(63) , x(64), x(67),x(68), x(71), x(75),x(78), x(79) , x(81), x(82) , x(83), x(84)...
                                            , x(86),x(88) , x(93),x(94), x(95) , x(96), x(97), x(98), x(101) , x(102) ,x(106), x(108), x(109) , x(114), x(115),x(117), x(118)]';% apply the optimal solution that was found 
                    %        Pd = Pd0 + delta_Pd_pu;
                             delta_Pd = delta_Pd_pu.*100; % Wasseem 07/18/2020, the base MVA is 100 , 

                             delta_Pg_pu = [x(1), x(4), x(6),	x(8),	x(10),	x(12),	x(15),	x(18),	x(19),	x(24), x(25),	x(26),	x(27),	x(31),	x(32),	x(34),	x(36),	x(40),	x(42),	x(46),	x(49) ...
                                            , x(54),	x(55),	x(56),	x(59),	x(61),	x(62),	x(65),	x(66),	x(69),	x(70),	x(72),	x(73),	x(74),	x(76),	x(77),	x(80),	x(85),	x(87),	x(89),	x(90),	x(91),	x(92) ...	
                                            , x(99),	x(100),	x(103),	x(104),	x(105),	x(107),	x(110),	x(111),	x(112),	x(113),	x(116)]'; 
                             delta_Pg = delta_Pg_pu.*100;

                    
                        else % this else is added to make sure that the function always outputs value !
                            % the following statemnt can be taken when the pointer
                     
                    %          delta_Pg = delta_Pg_pu .* 100;
                        end 

                    end

                    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
                    % build an optimization problem based on the sensitivity matrix 
                    %

                    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

                             delta_Pd_pu = [x(2),x(3) , x(5),x(7), x(9) , x(11), x(13), x(14), x(16) , x(17) ,x(20), x(21), x(22) , x(23), x(28),x(29), x(30), x(33),x(35), x(37) , x(38), x(39) , x(41), x(43)...
                                            , x(44),x(45) , x(47),x(48), x(50) , x(51), x(52), x(53) , x(57) ,x(58), x(60), x(63) , x(64), x(67),x(68), x(71), x(75),x(78), x(79) , x(81), x(82) , x(83), x(84)...
                                            , x(86),x(88) , x(93),x(94), x(95) , x(96), x(97), x(98), x(101) , x(102) ,x(106), x(108), x(109) , x(114), x(115),x(117), x(118)]'; % apply the optimal solution that was found 
                    %        Pd = Pd0 + delta_Pd_pu;
                             delta_Pd = delta_Pd_pu.*100; % Wasseem 07/18/2020, the base MVA is 100 , 
                             delta_Pg_pu = -1.*[x(1), x(4), x(6),	x(8),	x(10),	x(12),	x(15),	x(18),	x(19),	x(24), x(25),	x(26),	x(27),	x(31),	x(32),	x(34),	x(36),	x(40),	x(42),	x(46),	x(49) ...
                                            , x(54),	x(55),	x(56),	x(59),	x(61),	x(62),	x(65),	x(66),	x(69),	x(70),	x(72),	x(73),	x(74),	x(76),	x(77),	x(80),	x(85),	x(87),	x(89),	x(90),	x(91),	x(92) ...	
                                            , x(99),	x(100),	x(103),	x(104),	x(105),	x(107),	x(110),	x(111),	x(112),	x(113),	x(116)]'; 
                             delta_Pg = delta_Pg_pu.*100; 
                
                
                end
                
              
                    
                    
            
           
        
        
        
        end

end