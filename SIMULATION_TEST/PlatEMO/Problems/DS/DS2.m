classdef DS2 < PROBLEM
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
        function obj = DS2()
            K = 10;
            r = 0.25;
            alpha = 1;
            gamma = 4;
            tao = 1;
            obj.Parameter = table(K,r,alpha,gamma,tao);
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = K;
            obj.Global.l_D = K;
            
            obj.Global.u_lower    = [0.001,-K*ones(1,obj.Global.u_D-1)];
            obj.Global.u_upper    = K*ones(1,obj.Global.u_D);
            obj.Global.l_lower    = -K*ones(1,obj.Global.l_D);
            obj.Global.l_upper    = K*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
            K = obj.Parameter.('K');
            r = obj.Parameter.('r');
%             alpha = obj.Parameter.('alpha');
            gamma = obj.Parameter.('gamma');
            tao = obj.Parameter.('tao');
            
            N = size(uDecs,1);
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    v1=zeros(N,1);
                    v2=zeros(N,1);
                    y1 = y(:,1);
                    v1(y1>1) = y1(y1>1)-(1-cos(0.2*pi));
                    v1(~(y1>1)) = cos(0.2*pi)*y1(~(y1>1))+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    v2(y1>1) = 0.1*(y1(y1>1)-1)-sin(0.2*pi);
                    v2(~(y1>1)) = -sin(0.2*pi)*y1(~(y1>1))+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    
                    u_Obj(:,1) = v1+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    u_Obj(:,2) = v2+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                             
                    l_Obj = [];
                case 'lower'
                    u_Obj = [];
                    
                    l_Obj(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2);
                    l_Obj(:,2) = sum((x-y).^2.*repmat(1:K,N,1),2);
                case 'bilevel'
                    v1=zeros(N,1);
                    v2=zeros(N,1);
                    y1 = y(:,1);
                    v1(y1>1) = y1(y1>1)-(1-cos(0.2*pi));
                    v1(~(y1>1)) = cos(0.2*pi)*y1(~(y1>1))+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    v2(y1>1) = 0.1*(y1(y1>1)-1)-sin(0.2*pi);
                    v2(~(y1>1)) = -sin(0.2*pi)*y1(~(y1>1))+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    
                    u_Obj(:,1) = v1+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    u_Obj(:,2) = v2+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                             
                    l_Obj(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2);
                    l_Obj(:,2) = sum((x-y).^2.*repmat(1:K,N,1),2);         
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,parameters,type)
            N = size(uDecs,1);
            u_Con = zeros(N,1);
            l_Con = u_Con;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\DS2.mat']);
            catch                
                N = 1e4;
               
                r = obj.Parameter.('r');
                y1=[0.001 0.2 0.4 0.6 0.8 1];
                v1 = cos(0.2*pi)*y1+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
                v2 = -sin(0.2*pi)*y1+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
                t=(2:1/(N-1):3)';
                P(:,1) = reshape(cell2mat(cellfun(@(x) x+r*cos(pi/2*t), {v1}, 'UniformOutput',false)),[],1);
                P(:,2) = reshape(cell2mat(cellfun(@(x) x+r*sin(pi/2*t), {v2}, 'UniformOutput',false)),[],1);
                Front = NDSort(P,1);
                P = P(Front==1,:);
                
                if nargin >1
                [W,~] = UniformPoint(n,obj.Global.u_M);
                [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
                P = P(I,:);
                save([path,'\DS2.mat'],'P');
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