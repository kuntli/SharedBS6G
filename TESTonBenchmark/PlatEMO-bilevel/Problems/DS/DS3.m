classdef DS3 < PROBLEM
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
        function obj = DS3()
            K = 10;
            r = 0.2;
            tao = 1;
            obj.Parameter = table(K,r,tao);
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = K;
            obj.Global.l_D = K;
            
            obj.Global.u_lower    = zeros(1,obj.Global.u_D);
            obj.Global.u_upper    = K*ones(1,obj.Global.u_D);
            obj.Global.l_lower    = -K*ones(1,obj.Global.l_D);
            obj.Global.l_upper    = K*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
            K = obj.Parameter.('K');
%             r = obj.Parameter.('r');
            tao = obj.Parameter.('tao');
            
            N = size(uDecs,1);
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    y1 = y(:,1);
                    y2 = y(:,2);
                    R = 0.1+0.15*abs(sin(2*pi*(y1-0.1)));
                    temp = repmat(cell2mat(cellfun(@(x) x/2, {3:K}, 'UniformOutput',false)),N,1);
                    u_Obj(:,1) = y1+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*cos(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps)));
                    u_Obj(:,2) = y2+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*sin(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps)));   
                
                    l_Obj = [];
                case 'lower'
                    u_Obj = [];

                    l_Obj(:,1) = x(:,1) + sum((x(:,3:end)-y(:,3:end)).^2,2);
                    l_Obj(:,2) = x(:,2) + sum((x(:,3:end)-y(:,3:end)).^2,2);
                case 'bilevel'
                    y1 = y(:,1);
                    y2 = y(:,2);
                    R = 0.1+0.15*abs(sin(2*pi*(y1-0.1)));
                    temp = repmat(cell2mat(cellfun(@(x) x/2, {3:K}, 'UniformOutput',false)),N,1);
                    u_Obj(:,1) = y1+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*cos(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps)));
                    u_Obj(:,2) = y2+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*sin(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps))); 

                    l_Obj(:,1) = x(:,1) + sum((x(:,3:end)-y(:,3:end)).^2,2);
                    l_Obj(:,2) = x(:,2) + sum((x(:,3:end)-y(:,3:end)).^2,2);
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            switch type
                case 'upper'
                    u_Con =(1-uDecs(:,1).^2)-uDecs(:,2);

                    l_Con = [];
                case 'lower'
                    u_Con = [];

                    r = obj.Parameter.('r');                   
                    l_Con = (uDecs(:,1)-lDecs(:,1)).^2+(uDecs(:,2)-lDecs(:,2)).^2-r^2;
                case 'bilevel'
                    u_Con =(1-uDecs(:,1).^2)-uDecs(:,2);

                    r = obj.Parameter.('r');                   
                    l_Con = (uDecs(:,1)-lDecs(:,1)).^2+(uDecs(:,2)-lDecs(:,2)).^2-r^2;
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\DS3.mat']);
            catch

                N = 1e4;

                y1 = (0:0.1:1.3)';
                y2 = max(1-y1.^2,0);
                R = 0.1+0.15*abs(sin(2*pi*(y1-0.1)));
                t=1+0.5*(0:1/(N-1):1)';
                P(:,1) = cell2mat(cellfun(@(x,y) y+x.*cos(pi*t), num2cell(R), num2cell(y1),'UniformOutput',false));
                P(:,2) = cell2mat(cellfun(@(x,y) y+x.*sin(pi*t), num2cell(R), num2cell(y2),'UniformOutput',false));
                Front = NDSort(P,1);
                P = P(Front==1,:);
                
                if nargin >1
                [W,~] = UniformPoint(n,obj.Global.u_M);
                [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
                P = P(I,:);
                save([path,'\DS3.mat'],'P');
                end
            end
        end
        
        function P = lower_PF(obj,y)
            N = 10000;
            r = obj.Parameter.('r');
            t=(2:1/(N-1):3)';
            P(:,1) = y(:,1) + r*cos(pi/2*t);
            P(:,2) = y(:,2) + r*sin(pi/2*t);
        end
    end
end