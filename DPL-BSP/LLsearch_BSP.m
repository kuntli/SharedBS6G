function [SP,RP] = LLsearch_BSP(xu,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    Nl = Global.Nl;
    u_D = Global.u_D;
    l_D = Global.l_D;
    l_upper = Global.l_upper;
    l_lower = Global.l_lower;
    l_M = Global.l_M;
    u_M = Global.u_M;
    l_evaluated0 = Global.l_evaluated; 
    if Global.problem.Parameter.area<6000
        lFEs_max = 500;
    elseif Global.problem.Parameter.area<8000
        lFEs_max = 700;
    else
        lFEs_max = 1000;
    end
    
    Record = cell(1,10);
    Record_max = zeros(10,Global.l_M);  
    HV = zeros(1,10);
    
    if isa(xu(1),'INDIVIDUAL')
        Xls = round(rand(Nl-length(xu),l_D).*(l_upper-l_lower)+l_lower);
        Xls(rand(size(Xls))>0.5) = 0;
        Xls = reshape(Xls,size(Xls,1),u_D,l_D/u_D).*xu(1).udec;
        Xls = reshape(Xls,size(Xls,1),l_D);
        RP =[xu,INDIVIDUAL(repmat(xu(1).udec,size(Xls,1),1),Xls,'lower')];
        
    else
        Xl = round(rand(Nl,l_D).*(l_upper-l_lower)+l_lower);
        Xl(rand(size(Xl))>0.5) = 0;
        Xl = reshape(Xl,Nl,u_D,l_D/u_D).*xu; 
        Xl = reshape(Xl,Nl,l_D);
        RP = INDIVIDUAL(repmat(xu,Nl,1),Xl,'lower'); %ReferencePoints
    end
    Xu = RP.udecs;
   
    
    notermination = true;
    gen = 1;
    
    [Z,~] = UniformPoint(Nl,l_M);
    
    [W,~] = UniformPoint(Nl,u_M);
    
    [FrontNo,MaxFNo_sp] = NDSort(RP.lobjs,RP.lcons,Nl);
    [~,RP] = adds(RP,FrontNo');
%     CrowdDis = CrowdingDistance(RP.lobjs,FrontNo);
    
%     Next_SP = FrontNo < MaxFNo_sp;
%     Last     = find(FrontNo==MaxFNo_sp);
%     [~,Rank] = sort(CrowdDis(Last),'descend');
%     Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;
%     SP = RP(Next_SP);
%     FrontNo = FrontNo(Next_SP);
%     CrowdDis = CrowdDis(Next_SP);

    SP = RP;
    
    Record{1} = RP(RP.adds ==1);
    Record_max(1,:) = max(RP(RP.adds ==1).lobjs);
    ReferencePoint0 = max(Record_max);
    
%     l_PF = Global.problem.lower_PF(Xu(1,:));
    while notermination
        %% Display 
%         if gen>1
%             delete(subplot(2,2,2));
%             delete(subplot(2,2,4));
%         end
%         subplot(2,2,2);        
%         title('Lower-level Obj')
%         Draw(RP.lobjs);
% %         Draw(l_PF,'r');
%         I = sum(SP.lcons>0,2)<=0;
%         Draw(SP(I).lobjs,'go');
%         Draw(SP(~I).lobjs,'ro');        
%         
%         subplot(2,2,4);        
%         title('Lower-level Dec')
%         Draw(Xu(1,:),'b');
%         Draw(SP(I).ldecs,'g');
%         Draw(SP(~I).ldecs,'r');
%         
%         pause(0.0000001)
        %% Generate the offspring
%          MatingPool = TournamentSelection(2,Nl,FrontNo,-CrowdDis);
%          Pl = SP(MatingPool).ldecs;
         
        Pl = SP.ldecs;
        I = zeros(Nl,2);
        for i= 1:Nl
           I(i,:) =  randperm(Nl,2);
        end
        P1 = Pl(I(:,1),:);
        P2 = Pl(I(:,2),:);
        Ql = round(DE(Pl,P1,P2,'lower'));
    
%         I =randperm(Nl);
%         Ql = round(GA(Pl(I,:),'lower'));

        Ql = reshape(Ql,Nl,u_D,l_D/u_D).*xu; 
        Ql = reshape(Ql,Nl,l_D);
    
        SQ = INDIVIDUAL(Xu,Ql,'lower');        
       
        SR = [SP,SQ,RP(~ismember(RP.ldecs,SP.ldecs,'rows'))];
        
        [FrontNo,MaxFNo_rp] = NDSort(SR.lobjs,SR.lcons,Nl);
        [~,SR] = adds(SR,FrontNo');
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lobjs,FrontNo);
        
        
        %% Selection
        Next_RP = FrontNo < MaxFNo_rp;                         
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo_rp);
%         [~,Rank] = sort(CrowdDis(Last),'descend');
%         Next_RP(Last(Rank(1:Nl-sum(Next_RP)))) = true;

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
           
           
           
%             zmin = min(SR(all(SR.lcons<=0,2)).lobjs,[],1);
%             if isempty(zmin)
%                 zmin = min(SR.uobjs,[],1);
%             end
            
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
            
%             u_CrowdDis = CrowdingDistance(Candidate.uobjs,u_FrontNo);
%             last     = find(u_FrontNo==u_MaxFNo);
%             [~,rank] = sort(u_CrowdDis(last),'descend');
%             next(last(rank(1:Nl-sum(next)))) = true;
%             
            Next_SP = false(1,length(SR));
            Next_SP(index(next)) = true;           
       end
        RP = SR(Next_RP);
        SP = SR(Next_SP);
%         FrontNo = FrontNo(Next_SP);
%         CrowdDis = CrowdDis(Next_SP);
      
        %% termination check
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = RP(RP.adds==1);
        Record_max(k,:) = max(RP(RP.adds==1).lobjs);     
        
            ReferencePoint = max(Record_max);            
            
            if gen==2 || sum(ReferencePoint~=ReferencePoint0)>0
               ReferencePoint0 = ReferencePoint;               
               for i = 1 : min(gen,10)
                  HV(i) = hv(Record{i}.lobjs,ReferencePoint); 
               end
            else
                HV(k) = hv(Record{k}.lobjs,ReferencePoint);
            end
            
            HV_max = max(HV(1 : min(gen,10)));
            HV_min = min(HV(1 : min(gen,10)));      
        if gen >= 10              
            improverate = (HV_max - HV_min)/(HV_max + HV_min);
            tm = max(1,ceil(-log10(improverate))*10);
            notermination = improverate > 0.001;%&&Global.l_evaluated-l_evaluated0<lFEs_max;                     
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