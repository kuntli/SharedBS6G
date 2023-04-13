function [u_Decs,Archive,output] = Initilizing_ns(Global)
% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions
%
    %% Setting
    Nu = Global.Nu;
    Nl = Global.Nl;
    u_M = Global.u_M;
    
    
    %% Initilizing 
    u_Decs = Global.problem.Init(Nu,'upper');
    Archive = cell(1,Nu);
    output(1,Nu*Nl) = INDIVIDUAL;   
    Nr = 0;
    
    % Lower Level (LL) multiobjective optimization run for all xu,i
    for i = 1:Nu
        
        Archive{i}= LLsearch_ns(u_Decs(i,:),Global);           %Update SPi to LL Pareto-optimum

        output(Nr+1:Nr+sum(Archive{i}.adds==1)) = Archive{i}(Archive{i}.adds==1);
        Nr = Nr +sum(Archive{i}.adds==1);
  
    end
    output(Nr+1:end) = [];
    
              
      %% Update 
   [FrontNo,~] = NDSort(output.uobjs,[output.ucons,output.lcons],1);
   output = output(FrontNo==1);
   if length(output)>Nu*Nl
       [W,~] = UniformPoint(Nu*Nl,u_M);
       [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
       output = output(I);
   end
end