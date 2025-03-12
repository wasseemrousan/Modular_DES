

ps1= case300_002_ps;
ps2 = case300_2;
ps3 = case300; 
matrixmap(:,1) = [ps3.branch(:,1);ps3.branch(:,2)];
matrixmap(:,2) = [ps1.branch(:,1);ps1.branch(:,2)];

matrixmap2= unique(matrixmap(:,2),'stable');
matrixmap1= unique(matrixmap(:,1),'stable');
matrixmap_clean = [matrixmap1 matrixmap2];



psorig= ps2;
for e=1:300

idx1 = find(ps2.bus(:,1) == matrixmap_clean(e,1)); 

psorig.bus(idx1,1) = matrixmap_clean(e,2);

if ((psorig.bus(e,2) == 2) || (psorig.bus(e,2) == 3) )
    idx2 = find(ps2.gen(:,1) == ps2.bus(e,1)); 
    psorig.gen(idx2,1) = matrixmap_clean((find(ps2.gen(idx2,1) == matrixmap_clean(:,1))),2);
end

end 



save('psorig.mat','psorig');