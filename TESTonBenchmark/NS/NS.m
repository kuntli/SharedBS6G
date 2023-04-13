function NS(Global)
% <algorithm> <N>


%------------------------------- Reference --------------------------------
% W. Wang, H.-L. Liu, and H. Shi, ¡°A multi-objective bilevel optimisation
%evolutionary algorithm with dual populations lower-level search,¡±
%Connection Science, vol. 34, no. 1, pp. 1556¨C1581, 2022. [Online].
%Available: https://doi.org/10.1080/09540091.2022.2077312

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
     warning off
    
    
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.u_M);
    Hv = zeros(1,10);
    
    notermination = true;    
    gen = 1;
    

    u_PF = Global.problem.PF(1000);
    if ~isempty(u_PF)
        tRP = 1.1*max(u_PF);
    end
    
  tic;
    %% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions


    [u_Decs,Archive,output] = Initilizing_ns(Global);
   
    time = toc;
    times = time;
    Record{1} = output;
    Record_max(1,:) = max(output.uobjs);
    
    NFEs = [Global.u_evaluated,Global.l_evaluated];
    if isempty(u_PF)
        Igd  = []; 
        hv   = [];
    else
        Igd  = IGD(output.uobjs,u_PF);
        hv   = HV(output.uobjs,tRP);
    end
    outputs = cell(1,1);
    outputs{1} = output; 
    
    Output = cell(1,6);
    Output{1} = NFEs;
    Output{2} = output;
    Output{3} = Igd;
    Output{4} = hv;
    Output{5} = outputs;
    Output{6} = times;
    
    %% Optimization
    while Global.NotTermination(Output,notermination)
    tic;
        %% Generate the upper level Offspring
        SQ = Reproduce_ns(u_Decs,Global);
%% Combine parents and offspring, Update the Population
    
%         Draw(cat(2,Archive{:}).uobjs);
%         Draw(cat(2,SQ{:}).uobjs,'yo');
        [u_Decs,Archive,output] = UL_Select_ns(output,Archive,SQ,Global);   
%         Draw(output.uobjs,'go');    
%         pause(0.01);
        
        time = time + toc;
        times = [times;time];
        gen = gen + 1;
        %% Termination_check
        % Record the last 10 generations of populations
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = output;
        Record_max(k,:) = max(output.uobjs); 
         
        
        ReferencePoint = max(Record_max); 
        if sum(Hv)==0|| sum(ReferencePoint~=ReferencePoint0)>0
            ReferencePoint0 = ReferencePoint;
            for i = 1 : min(gen,10)
                Hv(i) = HV(Record{i}.uobjs,ReferencePoint); 
            end
         else
            Hv(k) = HV(Record{k}.uobjs,ReferencePoint);
         end
            
         HV_max = max(Hv(1 : min(gen,10)));
         HV_min = min(Hv(1 : min(gen,10)));
            

        if gen > 10    
            notermination = (HV_max - HV_min)/(HV_max + HV_min) > 0.001; 
        end
        
        NFEs = [NFEs;Global.u_evaluated,Global.l_evaluated];
        if ~isempty(u_PF)
            Igd  = [Igd;IGD(output.uobjs,u_PF)];
            hv   = [hv;HV(output.uobjs,tRP)];
        end
        newP = cell(1);
        newP{1} = output;
        outputs = [outputs,newP];
        
        
        Output{1} = NFEs;
        Output{2} = output;
        Output{3} = Igd;
        Output{4} = hv;
        Output{5} = outputs;
        Output{6} = times;
        
    end
end