function RP = LLsearch_ns(x,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    Nl = Global.Nl;
    Sn = [1:Nl];

    l_M = Global.l_M;
    u_M = Global.u_M;
    
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.l_M);  
    Hv = zeros(1,10);
    
    if isa(x(1),'INDIVIDUAL')
        if length(x)<Nl
        Xl = Global.problem.Init(Nl-length(x),'lower');
        RP =[x,INDIVIDUAL(repmat(x(1).udec,size(Xl,1),1),Xl,'lower')];
        else
            RP = x;
        end
    else       
        Xl = Global.problem.Init(Nl,'lower');
        RP = INDIVIDUAL(repmat(x,Nl,1),Xl,'lower'); %ReferencePoints
    end
    xu = RP(1).udecs;
   
    
    notermination = true;
    gen = 1;
    
    
    [FrontNo,~] = NDSort(RP.lobjs,RP.lcons,Nl);
    [~,RP] = adds(RP,FrontNo');
    CrowdDis = CrowdingDistance(RP.lobjs,FrontNo);
    
    Record{1} = RP(RP.adds ==1);
    Record_max(1,:) = max(RP(RP.adds ==1).lobjs);
    ReferencePoint0 = max(Record_max);
    
    while notermination
        %% Generate the offspring
        MatingPool = TournamentSelection(2,Nl,FrontNo,-CrowdDis);
        Pl = RP.ldecs;
        P1 = RP(MatingPool).ldecs;

        I = zeros(Nl,2);
        for i= 1:Nl
           Sn_temp = Sn;
           Sn_temp(i)=[];
           I(i,:) =  Sn_temp(randperm(Nl-1,2));
        end
        P2 = Pl(I(:,1),:);
        P3 = Pl(I(:,2),:);
        Ql = DE(P1,P2,P3,'lower');


        SQ = INDIVIDUAL(repmat(xu,size(Ql,1),1),Ql,'lower');        
       
        SR = [RP,SQ];
        
        [FrontNo,MaxFNo] = NDSort(SR.lobjs,SR.lcons,Nl);
        [~,SR] = adds(SR,FrontNo');
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lobjs,FrontNo);
        
        
        %% Selection
        Next_RP = FrontNo < MaxFNo;                         
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next_RP(Last(Rank(1:Nl-sum(Next_RP)))) = true;
        RP = SR(Next_RP);
        FrontNo = FrontNo(Next_RP);
        CrowdDis =CrowdDis(Next_RP);

        %% termination check
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = RP(RP.adds==1);
        Record_max(k,:) = max(RP(RP.adds==1).lobjs);     
        
            ReferencePoint = max(Record_max);            
            
            if sum(Hv)==0 || sum(ReferencePoint~=ReferencePoint0)>0
               ReferencePoint0 = ReferencePoint;               
               for i = 1 : min(gen,10)
                  Hv(i) = HV(Record{i}.lobjs,ReferencePoint); 
               end
            else
                Hv(k) = HV(Record{k}.lobjs,ReferencePoint);
            end
            
            HV_max = max(Hv(1 : min(gen,10)));
            HV_min = min(Hv(1 : min(gen,10)));      
        if gen >= 10              
            improverate = (HV_max - HV_min)/(HV_max + HV_min);
            
            notermination = improverate > 0.001;                     
        end
            
        gen = gen + 1;
            
            
    end
    
    index = find(RP.adds == 1);
    for i = 1:length(index)
        if isempty(RP(index(i)).uobj)
            RP(index(i)) = INDIVIDUAL(RP(index(i)),'upper');
        end
    end

end
