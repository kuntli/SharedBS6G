classdef TP1 < PROBLEM
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
        function obj = TP1()
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = 1;
            obj.Global.l_D = 2;
            
            obj.Global.u_lower    = zeros(1,obj.Global.u_D);
            obj.Global.u_upper    = ones(1,obj.Global.u_D);
            obj.Global.l_lower    = -ones(1,obj.Global.l_D);
            obj.Global.l_upper    = ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
           
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    u_Obj(:,1) =  x(:,1)-y;
                    u_Obj(:,2) = x(:,2);
                             
                    l_Obj = [];
                case 'lower'
                    u_Obj =[];
                    
                    l_Obj(:,1) = x(:,1);
                    l_Obj(:,2) = x(:,2);     
                case 'bilevel'
                    u_Obj(:,1) =  x(:,1)-y;
                    u_Obj(:,2) = x(:,2);
                             
                    l_Obj(:,1) = x(:,1);
                    l_Obj(:,2) = x(:,2);          
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            y = uDecs;
            x = lDecs;
            
            switch type
                case 'upper'
                    u_Con = -(1+x(:,1)+x(:,2));
                    l_Con = [];
                case 'lower'
                    u_Con = [];
                    l_Con = x(:,1).^2+ x(:,2).^2-y.^2;
                case 'bilevel'
                    u_Con = -(1+x(:,1)+x(:,2));
                    l_Con = x(:,1).^2+ x(:,2).^2-y.^2;
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            fullpath = mfilename('fullpath'); 
            [path,~]=fileparts(fullpath);
            try
                load([path,'\TP1.mat']);
            catch
                N = 1e6;
                y = (1:1/(N-1):sqrt(2))';
                x2 = [-1/2+1/4*sqrt(4*y.^2-4);-1/2-1/4*sqrt(4*y.^2-4)];
                x1 = -1-x2;
                P(:,1) = x1-1/sqrt(2)*[y;y];
                P(:,2) = x2;
                Z=min(P,[],1);
                W = UniformPoint(1e3,2);
                Next = false(1,size(P,1));
                for i=1:1e3
                    index = find(~Next);
                    [~,pos] = min(max((P(index,:)-Z)./repmat(W(i,:),length(index),1),[],2));
                    Next(index(pos)) = true;
                end
                P = P(Next,:);
                
                if nargin >1
                [W,~] = UniformPoint(n,obj.Global.u_M);
                [~,I] = max(1-pdist2(P-min(P),W,'cosine'),[],1);
                P = P(I,:);
                
                save([path,'\TP1.mat'],'P');
                end
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
            t=(2:1/(N-1):3)';
            P(:,1) = y*cos(pi/2*t);
            P(:,2) = y*sin(pi/2*t);
        end
    end
end