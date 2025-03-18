clc;
clear all;
C = psconstants; % tells me where to find my data

opt = psoptions;

ps = case300_002_ps;
ps = updateps(ps);
opt.sim.nHopLoc = 1; 
ps_agents = set_up_agents(ps,opt);




for t =1:300

  % we want to find the neighboring buses and all connected components
  ps_agents(t).loc_nei;
   idx = {}; col= {} ;  
  for e = 1:size(ps_agents(t).loc_nei,2)
      bid = ps_agents(t).loc_nei(e); % bus id of the neighborhood
  %+++++++++++++++start the opeartion for each bus in neighborhood++++++++++++++
          idx{e} = (ps_agents(bid).branch(:,[1 2]) ==  bid); 
          col{e} = find(idx{e});
  end      
            epty = []; 
            for w =1:size(col,2)
                epty = [ epty; cell2mat(col(w))];
                for m =1:size(epty,1)
                    if epty(m) > 516
                        epty(m) =  epty(m)-516;
                    end
                end
            end
            epty = unique(epty(:).')';
            if size(epty,1) > 13
                epty = col{1,(find(ps_agents(t).loc_nei==t))}; % take the lines attache directly to agent
                for m =1:size(epty,1)
                    if epty(m) > 516
                        epty(m) =  epty(m)-516;
                    end
                end

            end
            epty = unique(epty(:).')';

          for v= 1:size(epty,1)
              
        
              if epty(v) <= 9 
                    x{v} = ['line00' num2str(epty(v)) '.gen'];
        
              else 
                  if epty(v) <= 99 && epty(v) > 9
                    x{v} = ['line0' num2str(epty(v)) '.gen'];
        
                  else 
                    x{v} = ['line' num2str(epty(v)) '.gen'];
        
                  end
              end
          end  
        
        
          for n = 2:(size(epty,1))
               
              if  n == 2 
               paralleloperation300(x{1,n-1},x{1,n})
              else
               paralleloperation300(x{1,n},'ComComP.gen')
              end
              
          end
        
        
      for e = 1:size(ps_agents(t).loc_nei,2)
      % for e = 1:6
          % if e>9
          %       % break the for loop if e>11 
          %       break
          % end

          bid = ps_agents(t).loc_nei(e); % bus id
          if ~isempty(find(ps.gen(:,1)==bid, 1))
              
              idx2 = (ps.gen(:,1) == bid); 
              col2 = find(idx2);
              if col2 <= 9 
                    y = ['Genr00' num2str(col2) '.gen'];
        
              else 
                 if (col2 > 9) && (col2 <= 99)

                  
                    y = ['Genr0' num2str(col2) '.gen'];
                 else
                    y = ['Genr' num2str(col2) '.gen'];
                 end
              end
              paralleloperation300(y,'ComComP.gen'); % change name of 'ComCom' to something else
              
          end           
          if ~isempty(find(ps.shunt(:,1)==bid, 1))
              idx2 = (ps.shunt(:,1) == bid); 
              col2 = find(idx2);
              if col2 <= 9 
                    y = ['Load00' num2str(col2) '.gen'];
        
              else 
                  if (col2 > 9) && (col2 <= 99)
                    y = ['Load0' num2str(col2) '.gen'];
                  else
                    y = ['Load' num2str(col2) '.gen'];

                  end
              end
              paralleloperation300(y,'ComComP.gen'); % change name of 'ComCom' to something else
              
              
          end
            
      end
              if t <= 9 
                    newname{t} = ['Bus' '00' num2str(t) 'e.gen']; 
              else 
                if t > 9 && t <= 99

                    newname{t} = ['Bus' '0' num2str(t) 'e.gen'];
                else 
                    newname{t} = ['Bus' num2str(t) 'e.gen'];
                end
              end

              movefile('ComComP.gen',newname{t})

          clear x; 
          clear y;
 
end

