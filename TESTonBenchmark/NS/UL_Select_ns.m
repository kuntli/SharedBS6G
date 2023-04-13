function [u_Dec,Archive,output] = UL_Select_ns(output,Archive,SQ,Global)


    %% Setting
    Nu = Global.Nu;
    Nl = Global.Nl;
    u_Dec = zeros(Nu,Global.u_D);
    
    %% Loading
    SR = [Archive,SQ];
    R(1,length(SR)*Nl) = INDIVIDUAL;
    Idx_r = zeros(length(SR)*Nl,1);
    Nr = 0;
    for i = 1:length(SR)        
        R(Nr+1:Nr+sum(SR{i}.adds==1)) = SR{i}(SR{i}.adds==1);
        Idx_r(Nr+1:Nr+sum(SR{i}.adds==1)) = ones(sum(SR{i}.adds==1),1)*i;
        Nr = Nr +sum(SR{i}.adds==1);                   
    end
    R(Nr+1:end) = [];
    Idx_r(Nr+1:end) = [];
    
    
    %% Update 
    Q = cat(2,SQ{:});
    output = [output,Q(Q.adds == 1)];
    [FNo,~] = NDSort(output.uobjs,[output.ucons,output.lcons],1);
    output = output(FNo == 1);
    if length(output)>Nu*Nl
        [W,~] = UniformPoint(Nu*Nl,u_M);
        [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
        output = output(I);
    end
    
    
     [FNo,~] = NDSort(R.uobjs,[R.ucons,R.lcons],inf);
     CrowdDis = CrowdingDistance(R.uobjs,FNo);
     
    na = 0;
    maxFNo = 1;
    while na<Nu
     if sum(FNo==maxFNo)>0
         current = find(FNo==maxFNo);
         [~,I] = max(CrowdDis(current));
         u_Dec(na+1,:) = R(current(I)).udec;
         index = Idx_r(current(I));
         Archive{na+1} = SR{index};
         
         R(Idx_r==index) = [];
         FNo(Idx_r==index) = [];
         CrowdDis(Idx_r==index) = [];
         Idx_r(Idx_r==index) = [];
         na = na+1;
     else
         maxFNo = maxFNo +1;
     end
    end
end

