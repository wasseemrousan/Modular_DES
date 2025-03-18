

function [CptMatrix,al,sigmaFl] = Markovianlabelforcontrol(outputtreeObject,ps,pre_contingency_flows,br_outages)

% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% afunction to label Illegal states based on raltion number (33) in Markovian transition model paper
% We Will use the function contain to label thses states with 1's and get them back to the follwing steps in take control action function
% after this function is built, there will be no need for the labels at all
% steps 
% we  will also need to pass (ps) data into ths function to do the DC flow
% analysis 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
C = psconstants;
psMarkov = case118;
% psMarkov = updateps(psMarkov);
% psMarkov = rebalance(psMarkov);
% psMarkov = dcpf(psMarkov);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

% Patterns will be redifined as well 
% varible s here refers the lines states, 
% ps.branch(:,C.br.status) = 1;  
% ps.branch(br_outages,C.br.status) = 0;  

s = ps.branch(:,11)'; 

% lcrtitcal = argmin(al)
% al = ((FlMax) - (uFl))/ (sigmaFl)

% A: incidence matrix 
As = makeIncidence(psMarkov.bus, psMarkov.branch)';
% convert sparse matrix into full matrix 
A = full(As)'; 
%
% yl = 1/xl , 
y =  (1./ps.branch(:,4))'; % get the reactance of the lines, then find the admitance y for each line !
% update the adimatiance of lines from vector s 
y = (y.*s);
yt = diag(y); 
% preliminiry variables 
Atelda =  (sqrt(yt)*A); 
% B = pinv(A) returns the Moore-Penrose Pseudoinverse of matrix A.
% AteldaCross = pinv(Atelda); 
% Define means and covariance of injected power (from PS data)  
% uptcalculation = loadcase('case30');
% upt = [ uptcalculation.bus(:,3)./100];% WASSEEM : convert into per unit
% get the most recent means upt of the line flows 
% upt = [ ps.gen(:,2)./(100) ; -ps.shunt(:,2)./(100)];  % the ps modified case doesn't have 4 buses info so we rplaaced by zero (Wasseem 06/23/2020)

%% this part is modified to suite the 118 Bus system 
% upt = [psMarkov.bus(:,3)./100]; % this formula is used from the original definition of the case study from MATPOWER (to use arrangemnt of buses from 0 to 30)
upt = [(0.*psMarkov.bus(:,3))./100]; % this formula is modified for the IEEE 118 Bus system

for i=1:size(upt,1)
   if psMarkov.bus(i,2) == 1 % the type of bus 1 is a load bus
      upt(i) = (-psMarkov.bus(i,3))./100;
   end
end
% now assign values of power to the generators , this is modified for the
% IEEE 118 Bus system 
upt(10) = psMarkov.gen(5,2)/100;
upt(12) = (psMarkov.gen(6,2))/100; 
upt(25) = psMarkov.gen(11,2)/100;
upt(26) = psMarkov.gen(12,2)/100;
upt(31) = (psMarkov.gen(14,2))/100;
upt(46) = psMarkov.gen(20,2)/100;
upt(49) = psMarkov.gen(21,2)/100;
upt(54) = psMarkov.gen(22,2)/100;
upt(59) = psMarkov.gen(25,2)/100;
upt(61) = psMarkov.gen(26,2)/100;
upt(65) = psMarkov.gen(28,2)/100;
upt(66) = psMarkov.gen(29,2)/100;
upt(69) = psMarkov.gen(30,2)/100;
upt(80) = psMarkov.gen(37,2)/100;
upt(87) = psMarkov.gen(39,2)/100;
upt(89) = psMarkov.gen(40,2)/100;
upt(100) = psMarkov.gen(45,2)/100;
upt(103) = psMarkov.gen(46,2)/100;
upt(111) = psMarkov.gen(51,2)/100;



%% this part is for the covariance matrix, where the Covariance matrrix is con;  
% load('CptMatrix.mat')

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% we'll caculate the Cpt matrix here by making monte carlo simualtions
% using called function
% meanvaleu = ps.shunt(:,2)'; % assumin 
% deviation = 0.07.*meanvaleu; 
% the mean and variance will be retarcted from an excell  sheet inside the
% function 
CptMatrix = CalculateCovMatrix(ps); 
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

Cpt = CptMatrix.*(1/10000) ; % WASSEEM : convert into per unit
% We need to make sure of the operation T from the paper 
% Also make sure to make all the uFt entries positive, they did some
% procedure in the paper
% yt = yt';
global uFt;
uFt =  (sqrt( yt))*pinv(Atelda')*upt; 
CFt = (sqrt(yt))*pinv(Atelda')*Cpt*pinv(Atelda)*(sqrt(yt)); %  
% CFll is diagonal elememts of CFt Matrix ?? 
CFll = (diag(CFt)); 
% main relation dvariables 
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

% FlMax is the maximum rating of the transmission line 
% C = psconstants; % tells me where to find my data
% load case6ww_ps;
% psforFo = case6ww_ps ;
% psforFo.shunt(:,2) = ps.shunt(:,2); % this is to acomodate the variations in the loads, the 10 sampleas of the load 
% psforFo = dcpf(psforFo);
% 
% FlMaxo = psforFo.branch(:,C.br.Pf) % this is the rational flow distribution under normal operating condition
% FlMax = (FlMaxo./100)*1.2;
% this defintion of Flmax is from the original paper
% FlMax = (abs(pre_contingency_flows)./100).*1.2;
FlMax = ((ps.branch(:,8))./100);
% FlMax = abs(pre_contingency_flows./ps.branch(:,8));
% % post_contingency_flows = ps.branch(:,C.br.Pf);
% for i=1:size(pre_contingency_flows,1)
%     if (FlMax(i))< 2 
%         FlMax(i) = 2; 
%     end
% end
% % FlMax = FlMax' 

% FlMax = (ps.branch(:,8)./100);

sigmaFl = sqrt(CFll); 

% find mimum element of an array 
% I guess this element should be element by element devision 
al = ((FlMax) - abs(uFt))./ (sigmaFl); 



% save al;
%% After al is determined, we need to get the minimum of al, which means the shortest distance of a transmiison line to be tripped,  
[alcritical, Index] = min(al) ; %Wasseem 04/26/2020, lets ry and take the abs of al 
% we want to get an index for how much is the overload is for this minima %
% Wasseem 05/10/2020
overloaduantity = al(Index)*sigmaFl(Index);
% disp('al the whole array is :') 
% disp(al');
% disp('alcritical for control action is :') 
% disp (alcritical)
% disp('critical line for control action is :') 
% disp(Index)
% disp('overloaduantity for control action is :') 
% disp(overloaduantity)


% index will be the number of the line, which will be labled critical
% (illegal state in ouputtree object file) 

end