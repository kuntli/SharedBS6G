classdef DS1 < PROBLEM
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
        function obj = DS1()
            K = 10;
            r = 0.1;
            alpha = 1;
            gamma = 1;
            tao = 1;
            obj.Parameter = table(K,r,alpha,gamma,tao);
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = K;
            obj.Global.l_D = K;
            
            obj.Global.u_lower    = [1,-K*ones(1,obj.Global.u_D-1)];
            obj.Global.u_upper    = [4,K*ones(1,obj.Global.u_D-1)];
            obj.Global.l_lower    = -K*ones(1,obj.Global.l_D);
            obj.Global.l_upper    = K*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
            K = obj.Parameter.('K');
            r = obj.Parameter.('r');
            alpha = obj.Parameter.('alpha');
            gamma = obj.Parameter.('gamma');
            tao = obj.Parameter.('tao');
            
            N = size(uDecs,1);
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    temp = repmat(cell2mat(cellfun(@(x) (x-1)/2, {2:K}, 'UniformOutput',false)),N,1);
                    u_Obj(:,1) = (1+r-cos(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    u_Obj(:,2) = (1+r-sin(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                             
                    l_Obj = [];
                case 'lower'
                    u_Obj =[];
                    
                    l_Obj(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2)+10*sum(1-cos(pi/K*(x(:,2:end)-y(:,2:end))),2);
                    l_Obj(:,2) = sum((x-y).^2,2)+10*sum(abs(sin(pi/K*(x(:,2:end)-y(:,2:end)))),2);     
                case 'bilevel'
                    temp = repmat(cell2mat(cellfun(@(x) (x-1)/2, {2:K}, 'UniformOutput',false)),N,1);
                    u_Obj(:,1) = (1+r-cos(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    u_Obj(:,2) = (1+r-sin(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                             
                    l_Obj(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2)+10*sum(1-cos(pi/K*(x(:,2:end)-y(:,2:end))),2);
                    l_Obj(:,2) = sum((x-y).^2,2)+10*sum(abs(sin(pi/K*(x(:,2:end)-y(:,2:end)))),2);          
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            N = size(uDecs,1);
            u_Con = zeros(N,1);
            l_Con = u_Con;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\DS1.mat']);
            catch
                N = 1e6;
                r = obj.Parameter.('r');
                t=(2:1/(N-1):3)';
                P(:,1) = 1+r+(1+r)*cos(pi/2*t);
                P(:,2) = 1+r+(1+r)*sin(pi/2*t);
                if nargin >1
                [W,~] = UniformPoint(n,obj.Global.u_M);
                [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
                P = P(I,:);
                save([path,'\DS1.mat'],'P');
                end
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
%             r = obj.Parameter.('r');
            t=(0:1/(N-1):y(:,1))';
            P(:,1) = t.^2;
            P(:,2) = (t-y(:,1)).^2;
        end
    end
end