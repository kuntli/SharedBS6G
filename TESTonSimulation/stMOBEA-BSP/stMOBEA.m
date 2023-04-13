function stMOBEA(Global)
% <algorithm> <N>
% multiobjective bilevel evolutionary algorithm based on multiple quadratic fibers (mf-BLEA_CGAN)

%------------------------------- Reference --------------------------------

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
    
    Nu = Global.Nu;
    Nl = Global.Nl;
    K = 5;
    
    if Global.problem.Parameter.area<6000
        uFEs_max = 2500;
    elseif Global.problem.Parameter.area<8000
        uFEs_max = 500;
    else
        uFEs_max = 700;
    end
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.u_M);
    HV = zeros(1,10);
    
    notermination = true;    
    gen = 1;
    
    N0 = 5;
    N = N0;

    tmax = 40;
    tm = 1;
    %% Initialize upper level (UL) population and search for their optimal lower level (LL) solutions


    [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,W] = Initilizing_st(K,Global);
    TrainSet = l_Pop;
   
    Record{1} = output;
    Record_max(1,:) = max(output.uobjs);
    
    NFEs = [Global.u_evaluated,Global.l_evaluated];
    Igd  = []; %IGD(output.uobjs,u_PF);
    Hv   = [];%hv(output.uobjs,tRP);
    outputs = cell(1,1);
    outputs{1} = output; 
    
    Output = cell(1,5);
    Output{1} = NFEs;
    Output{2} = output;
    Output{3} = Igd;
    Output{4} = Hv;
    Output{5} = outputs;
    %% Optimization
    while Global.NotTermination(Output,notermination)
        %% train model
        planemodel = Plane_model(TrainSet);
        

        %% Generate the upper level Offspring
        
        [SQ,RQs] = Reproduce_st(tmax,u_Pop,Idx_ul,l_Pop,Arch_partion,TrainSet,planemodel,W,Global);
        TrainSet = [l_Pop,SQ];
        
        %% Combine parents and offspring, Update the Population
        if ~isempty(SQ)                          
            [u_Pop,Idx_ul,l_Pop,Arch_partion,output,Idx_oA,Archive,RPs,~] = UL_Select_BSP(l_Pop,output,Archive,RPs,SQ,RQs,W,Global);
        end    
            
        %% Refinesearch
        
        [output,Idx_oA,Archive,RPs] = RefineSearch_st(output,Idx_oA,Archive,RPs,tm,Global);

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
        if sum(HV)==0|| sum(ReferencePoint~=ReferencePoint0)>0
            ReferencePoint0 = ReferencePoint;
            for i = 1 : min(gen,10)
                HV(i) = hv(Record{i}.uobjs,ReferencePoint); 
            end
         else
            HV(k) = hv(Record{k}.uobjs,ReferencePoint);
         end
            
         HV_max = max(HV(1 : min(gen,10)));
         HV_min = min(HV(1 : min(gen,10)));
            
%          tm = max(1,ceil(-log10(max((HV_max - HV_min)/(HV_max + HV_min),0.001)))*15);
%          if gen >2
            index =mod([k+9,k+8],10);
            site  = index ==0;
            index(site) =10;
            
            Ir1 = max((HV(k)-HV(index(1)))/(HV_max + HV_min),0);
            Ir2 = max((HV(index(1))-HV(index(2)))/(HV_max + HV_min),0);
            if Ir2 >0.001&&Ir1>0.001
                tmax = min(40,max(floor(tmax + tmax/5*min(5,log2(Ir1/Ir2))),1));
            else
                tmax = 1;
            end
%          end
        if gen > 10    
            notermination = (HV_max - HV_min)/(HV_max + HV_min) > 0.001&&gen<=20;% && Global.u_evaluated<uFEs_max; 
            N =max(N,N0 + max(0,ceil(-log10(max((HV_max - HV_min)/(HV_max + HV_min),0.001)))*5));
        end
        
        NFEs = [NFEs;Global.u_evaluated,Global.l_evaluated];
%         Igd  = [Igd;IGD(output.uobjs,u_PF)];
%         Hv   = [Hv;hv(output.uobjs,tRP)];
        newP = cell(1);
        newP{1} = output;
        outputs = [outputs,newP];
        
        
        Output{1} = NFEs;
        Output{2} = output;
        Output{3} = Igd;
        Output{4} = Hv;
        Output{5} = outputs;
        %% Display
        
%         if gen>2
%             delete(subplot(2,2,1));
%             delete(subplot(2,2,3));
%         end
%         subplot(2,2,1);       
%         title('Upper-level Obj')
%         for i = 1:Nu
%             Draw(l_Pop{i}.uobjs);
%         end
% %         Draw(u_PF,'rs');
%         I = sum(u_Pop.ucons,2)<=0;
%         Draw(u_Pop(I).uobjs,'go');
%         Draw(u_Pop(~I).uobjs,'ro');
%         
%         
%         subplot(2,2,3);        
%         title('Upper-level Dec')
%         for i = 1:Nu
%             Draw(l_Pop{i}(1).udec);
%         end
%         if Global.u_D>3
%             for i = 1:Nu
%                 Draw(l_Pop{i}(1).udec);
%             end
%             Draw(u_Pop(I).udecs,'g');
%             Draw(u_Pop(~I).udecs,'r');
%         else
%             for i = 1:Nu
%                 plot(l_Pop{i}(1).udec,'o');
%             end
%             plot(u_Pop(I).udecs,'go');
%             plot(u_Pop(~I).udecs,'ro');
%         end
%         hold off;
%         pause(0.00001);
        
    end
end