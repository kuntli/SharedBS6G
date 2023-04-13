function SQ = Reproduce_ns(u_Decs,Global)

%% setting
    Nu = Global.Nu;
    Nl = Global.Nl;
    Sn = [1:Nu];
    
    Idx = zeros(Nu,2);
    for i = 1:Nu
        Sn_temp = Sn;
        Sn_temp(i)=[];
        Idx(i,:) = Sn_temp(randperm(Nu-1,2));
    end
    P2 = u_Decs(Idx(:,1),:);
    P3 = u_Decs(Idx(:,2),:);
     
    Qu = DE(u_Decs,P2,P3,'upper');
                             
    I = ismember(Qu,u_Decs,'rows');
    Qu(I,:) =[];
    Qu = unique(Qu,'rows');
    
    SQ = cell(1,size(Qu,1));
    for i = 1:size(Qu,1)
        SQ{i}= LLsearch_ns(Qu(i,:),Global); 
    end
end

