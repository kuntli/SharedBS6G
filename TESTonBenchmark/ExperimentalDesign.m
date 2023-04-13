function ExperimentalDesign
%% Settings
clc
clear
close all
%获取路径
fullpath = mfilename('fullpath'); 
[path,~]=fileparts(fullpath);

Need_Experiment   =   0;
Need_Statistics   =   1;
Need_Drawing      =   1;
Need_timecoputing =   0;
 
%运行次数
Run_Times     =    21;
%测试问题
Problems = {...   
                'TP1',...
                'TP2',...
                'TP3',...
                'DS1',...
                'DS2',...
                'DS3',...
                'DS4',...
                'DS5',...  
                'TP4'
           };
%测试算法
Algorithms = {...

            'stMOBEA'
            'BLEA_PM'
            'SABLEA'
            'MOBEA_DPL'
            'NS'
    };
Algorithm_name = {...
            'SABLEA-PM'
            'BLEA-PM'
            'SABLEA'           
            'MOBEA-DPL'
            'NS'
    };
 
  
%上层目标维数
u_M        = [...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2
    ];
%下层目标维数
l_M        = [...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2,...
                2
    ];
%上层自变量空间维数
u_D        = [...
                1,...
                1,...
                1,...
                10,...
                10,...
                10,...
                1,...
                1,...
                2
    ];
%下层自变量空间维数
l_D        = [...
                2,...
                14,...
                2,...
                10,...
                10,...
                10,...
                9,...
                9,...
                3
    ];


%HV中位数的序号 
Idx     =[...
          11,15,19,20,11,21,10,0;...
          15,18, 8, 3, 7, 5, 3,0;...
           7, 3,17,10,21,13, 4,0;...
           9, 9,17, 7,20,17, 7,0;...
           6, 6, 5,10, 6,10,14,0;...
           6, 9, 7, 1, 1,20, 2,0;...  
];
 
%% Experiment

Pool = {};
if Need_Experiment
    for a=1:length(Algorithms)
        TheAlgorithm = eval(['@',Algorithms{a}]);
        for p = 1:1:length(Problems)
            TheProblem = eval(['@',Problems{p}]);
            for j=1:Run_Times
               filenanme = [Algorithms{a},'_',Problems{p},'_uM',num2str(u_M(p)),'_lM',num2str(l_M(p)),'_uD',num2str(u_D(p)),'_lD',num2str(l_D(p)),'_',num2str(j),'.mat'];
%                 if ~exist(filenanme)
                    Inputs = {'-algorithm',TheAlgorithm,'-problem',TheProblem,'-Nu',20,'-Nl',20,'-run',j,'-save',1};
                    Pool = cat(1,Pool,Inputs);
%                 end
            end
        end
    end
end

%并行开关
if ~isempty(Pool)
    Isparallel = 1;
    parallelnum = 4;
    if ~Isparallel
        for k=1:size(Pool,1)
            main(Pool{k,:},folder);
        end
    else
        parpool('local',parallelnum)
        parfor k=1:size(Pool,1)
            Notcompelted = true; 
            while Notcompelted
                try
                    main(Pool{k,:});
                    Notcompelted =false;
                catch
                    
                    
                end
            end
        end
        delete(gcp);
    end
end

%% 统计实验结果
referencepoints = zeros(length(Problems),2);
for p = 1:length(Problems)    
    try
        load([Problems{p},'.mat']);
        uPF = P;
        referencepoints(p,:) = 1.1*max(uPF-min(uPF))+min(uPF);
    catch
         maxObj_a = zeros(length(Algorithms),u_M(p));
            for a=1:length(Algorithms)
                maxObj_r = zeros(Run_Times,u_M(p));
                for r = 1:Run_Times
                   filenanme = [Algorithms{a},'_',Problems{p},'_uM',num2str(u_M(p)),'_lM',num2str(l_M(p)),'_uD',num2str(u_D(p)),'_lD',num2str(l_D(p)),'_',num2str(r),'.mat'];

                    load(filenanme);
                    Pop =  result{3}{2};
                    %Draw(Pop.uobjs,symbol2{a});
                    maxObj_r(r,:) = max(Pop.uobjs); 
                end
                maxObj_a(a,:) = max(maxObj_r);
            end
            referencepoints(p,:) = max(maxObj_a);
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Need_Statistics
    for a=1:length(Algorithms)
        X = [];
        for p = 1:length(Problems)
            uFEs =zeros(Run_Times,1);
            lFEs =zeros(Run_Times,1);
            IGD = zeros(Run_Times,1);
            Hv = zeros(Run_Times,1);
            
            try
                load([Problems{p},'.mat']);
                uPF = P;
            catch
                uPF =[];
            end
            referencepoint = referencepoints(p,:);
            
            for r = 1:Run_Times
                filenanme = [Algorithms{a},'_',Problems{p},'_uM',num2str(u_M(p)),'_lM',num2str(l_M(p)),'_uD',num2str(u_D(p)),'_lD',num2str(l_D(p)),'_',num2str(r),'.mat'];
                load(filenanme);
                uFEs(r) = result{1};
                lFEs(r) = result{2};
                
                Pop =  result{3}{2};
                if ~isempty(uPF)
                    IGD(r) = result{3}{3}(end);                   
                else
                    IGD(r) = nan;
                end
                Hv(r)= hv(Pop.uobjs,referencepoint);
            end
            X =[X,uFEs,lFEs,IGD,Hv];
            HV_med = median(Hv);
            Ic = find(Hv==HV_med);
            if length(Ic)>1
                I = randi(length(find(Hv==HV_med)));
                Idx(a,p) = Ic(I);
            else
                Idx(a,p) = Ic;
            end
        end
%         xlswrite([path,'\Results\results.xls'],X,Algorithm_name{a},'B3')
         xlswrite([path,'\Results\NS.xls'],X,Algorithm_name{a},'B3')
    end
end


%%  作图
if Need_Drawing
    symbol ={...
        'ro-',...
        'k*-',...
        'bs-',...
        'cd-',...
        'g^-',...
        'mx-',...
        };
    
    for p = 7:8%1:length(Problems)
        if p ==3||p==9
           continue; 
        end
        referencepoint = referencepoints(p,:);
        hold on
        for a=1:length(Algorithms)-1
            filenanme = [Algorithms{a},'_',Problems{p},'_uM',num2str(u_M(p)),'_lM',num2str(l_M(p)),'_uD',num2str(u_D(p)),'_lD',num2str(l_D(p)),'_',num2str(Idx(a,p)),'.mat'];
            load(filenanme);
            
            FEs =[0; sum(result{3}{1},2)];
            Hv = zeros(length(FEs),1);
%             FEs = sum(result{3}{1},2);
%             Hv = result{3}{4};
            for i = 1:length(FEs)-1
                if p == 3|| p==9
                    Hv(i+1) = HV(result{3}{5}{i}.uobjs,referencepoint);
                else
                    Hv(i+1) = hv(result{3}{5}{i}.uobjs,referencepoint);
                end
            end
            %            Hv =log(Hv);
            if (p == 7|| p==8)&&a ==5
%                 Hv = log(Hv);
                site =1:3;
            else
            site =1:5:length(Hv);
            end
            plot(FEs(site),Hv(site),symbol{a});
        end
        legend(Algorithm_name);
        xlabel('FEs');
%         if p == 7|| p==8
%             ylabel('log(HV)');
%         else
        ylabel('HV');
%         end
        savefig( [path,'\Results\Fig\',Problems{p}, '.fig']);
        clf;
        hold off
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Need_timecoputing
    for a=1:length(Algorithms)
        X = [];
        for p = 1:length(Problems)
            T =zeros(Run_Times,1);
            for r = 1:Run_Times
                filenanme = [Algorithms{a},'_',Problems{p},'_uM',num2str(u_M(p)),'_lM',num2str(l_M(p)),'_uD',num2str(u_D(p)),'_lD',num2str(l_D(p)),'_',num2str(r),'.mat'];
                load(filenanme);
                T(r) = result{4};
            end
            X = [X,T];
        end
        xlswrite([path,'\Results\times.xls'],X,Algorithm_name{a},'B3')
    end
end
end


