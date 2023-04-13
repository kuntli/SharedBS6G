function [SP,RP] = LLsearch_DPL(xu,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    Nl = Global.Nl;

    l_M = Global.l_M;
    u_M = Global.u_M;
    
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.l_M);  
    Hv = zeros(1,10);
    
    if isa(xu(1),'INDIVIDUAL')
        if length(xu)<Nl
        Xl = Global.problem.Init(Nl-length(xu),'lower');
        RP =[xu,INDIVIDUAL(repmat(xu(1).udec,size(Xl,1),1),Xl,'lower')];
        else
            RP = xu;
        end
    else       
        Xl = Global.problem.Init(Nl,'lower');
        RP = INDIVIDUAL(repmat(xu,Nl,1),Xl,'lower'); %ReferencePoints
    end
    xu = RP(1).udecs;
   
    
    notermination = true;
    gen = 1;
    
    [Z,~] = UniformPoint(Nl,l_M);
    
    [W,~] = UniformPoint(Nl,u_M);
    
    [FrontNo,MaxFNo_sp] = NDSort(RP.lobjs,RP.lcons,Nl);
    [~,RP] = adds(RP,FrontNo');

    SP = RP;
    
    Record{1} = RP(RP.adds ==1);
    Record_max(1,:) = max(RP(RP.adds ==1).lobjs);
    ReferencePoint0 = max(Record_max);
    
    while notermination
        %% Generate the offspring
         
        Pl = SP.ldecs;
        I = zeros(Nl,2);
        for i= 1:Nl
           I(i,:) =  randperm(Nl,2);
        end
        P1 = Pl(I(:,1),:);
        P2 = Pl(I(:,2),:);
        Ql = DE(Pl,P1,P2,'lower');


        SQ = INDIVIDUAL(repmat(xu,size(Ql,1),1),Ql,'lower');        
       
        SR = [SP,SQ,RP(~ismember(RP.ldecs,SP.ldecs,'rows'))];
        
        [FrontNo,MaxFNo_rp] = NDSort(SR.lobjs,SR.lcons,Nl);
        [~,SR] = adds(SR,FrontNo');
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lobjs,FrontNo);
        
        
        %% Selection
        Next_RP = FrontNo < MaxFNo_rp;                         
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo_rp);


        Zmin = min(SR(all(SR.lcons<=0,2)).lobjs,[],1);
        if isempty(Zmin)
            Zmin = min(SR.lobjs,[],1);
        end
        Choose = LastSelection(SR(Next_RP).lobjs,SR(Last).lobjs,Nl-sum(Next_RP),Z,Zmin);
        Next_RP(Last(Choose)) = true;
        
        

       MaxFNo_sp = 1;
       while sum(FrontNo<=MaxFNo_sp)<Nl
           MaxFNo_sp = MaxFNo_sp+1;
       end
        
        
       if MaxFNo_sp>1 ||sum(FrontNo <= 1)== Nl           
           Next_SP = FrontNo < MaxFNo_sp;
           Last     = find(FrontNo==MaxFNo_sp);
           [~,Rank] = sort(CrowdDis(Last),'descend');
           Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;
                     
            
       else
            index =find(FrontNo==1);
            for i = 1:length(index)
               if isempty(SR(index(i)).uobj)
                   SR(index(i)) = INDIVIDUAL(SR(index(i)),'upper');
               end
            end
            Candidate = SR(index);
            
            [u_FrontNo,u_MaxFNo] = NDSort(Candidate.uobjs,[Candidate.ucons,Candidate.lcons],Nl);
            next = u_FrontNo<u_MaxFNo;
            
            Zmin = min(Candidate(all(Candidate.ucons<=0,2)).uobjs,[],1);
            if isempty(Zmin)
                Zmin = min(Candidate.uobjs,[],1);
            end
            last     = find(u_FrontNo==u_MaxFNo);
            Choose = LastSelection(Candidate(next).uobjs,Candidate(last).uobjs,Nl-sum(next),W,Zmin);
            next(last(Choose)) = true;
                    
            Next_SP = false(1,length(SR));
            Next_SP(index(next)) = true;           
       end
        RP = SR(Next_RP);
        SP = SR(Next_SP);
      
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
    
    index = find(SP.adds == 1);
    for i = 1:length(index)
        if isempty(SP(index(i)).uobj)
            SP(index(i)) = INDIVIDUAL(SP(index(i)),'upper');
        end
    end
    
    [site,Loc] = ismember(RP.ldecs,SP.ldecs,'rows');
    RP(site) = SP(Loc(site));
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end