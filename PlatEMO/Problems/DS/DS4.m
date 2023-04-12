classdef DS4 < PROBLEM
% <problem> <DS>
% Benchmark multi-objective bi-level problem

%------------------------------- Reference --------------------------------
% Deb K , Sinha A . An Evolutionary Approach for Bilevel Multi-objective 
% Problems[J]. Communications in Computer & Information Science, 2009,
% 35:17-24.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function obj = DS4()
            K = 5;
            L = 4;

            obj.Parameter = table(K,L);
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = 1;
            obj.Global.l_D = K+L;
            
            obj.Global.u_lower    = ones(1,obj.Global.u_D);
            obj.Global.u_upper    = 2*ones(1,obj.Global.u_D);
            obj.Global.l_lower    = [0,-(K+L)*ones(1,obj.Global.l_D-1)];
            obj.Global.l_upper    = [1,(K+L)*ones(1,obj.Global.l_D-1)];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
            K = obj.Parameter.('K');
            L = obj.Parameter.('L');
            
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    y1 = y(:,1);
                    x1 = x(:,1);
                    temp = sum(x(:,2:K).^2,2)+1;
                    
                    u_Obj(:,1) = (1-x1).*temp.*y1;
                    u_Obj(:,2) = x1.*temp.*y1;

                    l_Obj = [];
                case 'lower'
                    u_Obj = [];

                    y1 = y(:,1);
                    x1 = x(:,1);
                    temp = sum(x(:,K+1:K+L).^2,2)+1;
                    l_Obj(:,1) = (1-x1).*temp.*y1;
                    l_Obj(:,2) = x1.*temp.*y1;
                case 'bilevel'
                    y1 = y(:,1);
                    x1 = x(:,1);
                    temp = sum(x(:,2:K).^2,2)+1;
                    
                    u_Obj(:,1) = (1-x1).*temp.*y1;
                    u_Obj(:,2) = x1.*temp.*y1;
                    
                    temp = sum(x(:,K+1:K+L).^2,2)+1;
                    l_Obj(:,1) = (1-x1).*temp.*y1;
                    l_Obj(:,2) = x1.*temp.*y1;                
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            switch type
                case 'upper'
                    u_Con = -((1-lDecs(:,1)).*uDecs(:,1)+1/2*lDecs(:,1).*uDecs(:,1)-1);

                    l_Con = [];
                case 'lower'
                    u_Con = [];

                    N = size(uDecs,1);                  
                    l_Con = zeros(N,1);
                case 'bilevel'
                    u_Con = -((1-lDecs(:,1)).*uDecs(:,1)+1/2*lDecs(:,1).*uDecs(:,1)-1);
                    N = size(uDecs,1);                  
                    l_Con = zeros(N,1);
            end
        end
        %% Sample reference points on Pareto front
         function P = PF(obj,n)
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\DS4.mat']);
            catch    
                N = 1e6;
 
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = -2*P(:,1)+2;
                
                if nargin >1
                [W,~] = UniformPoint(n,obj.Global.u_M);
                [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
                P = P(I,:);
                
                save([path,'\DS4.mat'],'P');
                end
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
            y1 = y(:,1);
            t=(0:1/(N-1):y1)';
            P(:,1) = t;
            P(:,2) = y1-t;
        end
    end
end