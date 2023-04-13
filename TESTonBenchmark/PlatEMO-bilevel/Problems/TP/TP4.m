classdef TP4 < PROBLEM
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
        function obj = TP4()
           
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = 2;
            obj.Global.l_D = 3;
            
            obj.Global.u_lower    = zeros(1,obj.Global.u_D);
            obj.Global.u_upper    = 1.5e2*ones(1,obj.Global.u_D);
            obj.Global.l_lower    = zeros(1,obj.Global.l_D);
            obj.Global.l_upper    = 1.5e2*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
            
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    u_Obj(:,1) = -(y(:,1)+9*y(:,2)+ 10*x(:,1)+x(:,2)+3*x(:,3));
                    u_Obj(:,2) = -(9*y(:,1)+2*y(:,2)+ 2*x(:,1)+7*x(:,2)+4*x(:,3));
                    
                    l_Obj =[];
                case 'lower'
                    u_Obj =[];
                    
                    l_Obj(:,1) = -(4*y(:,1)+6*y(:,2)+ 7*x(:,1)+4*x(:,2)+8*x(:,3));
                    l_Obj(:,2) = -(6*y(:,1)+4*y(:,2)+ 8*x(:,1)+7*x(:,2)+4*x(:,3));                    
                case 'bilevel'
                    u_Obj(:,1) = -(y(:,1)+9*y(:,2)+ 10*x(:,1)+x(:,2)+3*x(:,3));
                    u_Obj(:,2) = -(9*y(:,1)+2*y(:,2)+ 2*x(:,1)+7*x(:,2)+4*x(:,3));
                    
                    l_Obj(:,1) = -(4*y(:,1)+6*y(:,2)+ 7*x(:,1)+4*x(:,2)+8*x(:,3));
                    l_Obj(:,2) = -(6*y(:,1)+4*y(:,2)+ 8*x(:,1)+7*x(:,2)+4*x(:,3));
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            y = uDecs;
            x = lDecs;

            switch type
                case 'upper'
                    u_Con(:,1) = 3*y(:,1)+9*y(:,2)+ 9*x(:,1)+5*x(:,2)+3*x(:,3)-1039;
                    u_Con(:,2) = -4*y(:,1)-y(:,2)+ 3*x(:,1)-3*x(:,2)+2*x(:,3)-94;
                    
                    l_Con = [];
                case 'lower'
                    u_Con = [];
                    
                    l_Con(:,1) = 3*y(:,1)-9*y(:,2)-9*x(:,1)-4*x(:,2)-61;
                    l_Con(:,2) = 5*y(:,1)+9*y(:,2)+ 10*x(:,1)-x(:,2)-2*x(:,3)-924;
                    l_Con(:,3) = 3*y(:,1)-3*y(:,2)+x(:,2)+5*x(:,3)-420;
                case 'bilevel'
                    u_Con(:,1) = 3*y(:,1)+9*y(:,2)+ 9*x(:,1)+5*x(:,2)+3*x(:,3)-1039;
                    u_Con(:,2) = -4*y(:,1)-y(:,2)+ 3*x(:,1)-3*x(:,2)+2*x(:,3)-94;
                    
                    l_Con(:,1) = 3*y(:,1)-9*y(:,2)-9*x(:,1)-4*x(:,2)-61;
                    l_Con(:,2) = 5*y(:,1)+9*y(:,2)+ 10*x(:,1)-x(:,2)-2*x(:,3)-924;
                    l_Con(:,3) = 3*y(:,1)-3*y(:,2)+x(:,2)+5*x(:,3)-420;
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            y = [146.2955,28.9394];
            x = [0,67.9318,0];
            P(:,1) = -(y(:,1)+9*y(:,2)+ 10*x(:,1)+x(:,2)+3*x(:,3));
            P(:,2) = -(9*y(:,1)+2*y(:,2)+ 2*x(:,1)+7*x(:,2)+4*x(:,3));
        end
        
        function P = lower_PF(obj,y)
            P = [];
        end
    end
end