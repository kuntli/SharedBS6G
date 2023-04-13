classdef TP3 < PROBLEM
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
        function obj = TP3()
            
            obj.Global.u_M = 2;
            obj.Global.l_M = 2;
            obj.Global.u_D = 1;
            obj.Global.l_D = 2;
            
            obj.Global.u_lower    = zeros(1,obj.Global.u_D);
            obj.Global.u_upper    = 10*ones(1,obj.Global.u_D);
            obj.Global.l_lower    = zeros(1,obj.Global.l_D);
            obj.Global.l_upper    = 10*ones(1,obj.Global.l_D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function [u_Obj,l_Obj,paratmeters] = CalObj(obj,uDecs,lDecs,paratmeters,type)
           
            y = uDecs;
            x = lDecs;
            switch type
                case 'upper'
                    u_Obj(:,1) = x(:,1)+x(:,2).^2 + y(:,1) + sin(x(:,1)+y(:,1)).^2;
                    u_Obj(:,2) = cos(x(:,2)).*(0.1+y(:,1)).*exp(-x(:,1)./(0.1+x(:,2)));
                             
                    l_Obj = [];
                case 'lower'
                    u_Obj =[];
                    
                    l_Obj(:,1) = ((x(:,1)-2).^2+(x(:,2)-2).^2)/4 + (x(:,2).*y(:,1)+(5-y(:,1)).^2)/16+sin(x(:,2)/10);
                    l_Obj(:,2) = (x(:,1).^2+(x(:,2)-6).^4-2*x(:,1).*y(:,1)-(5-y(:,1)).^2)/80;
                case 'bilevel'
                    u_Obj(:,1) = x(:,1)+x(:,2).^2 + y(:,1) + sin(x(:,1)+y(:,1)).^2;
                    u_Obj(:,2) = cos(x(:,2)).*(0.1+y(:,1)).*exp(-x(:,1)./(0.1+x(:,2)));
                             
                    l_Obj(:,1) = ((x(:,1)-2).^2+(x(:,2)-2).^2)/4 + (x(:,2).*y(:,1)+(5-y(:,1)).^2)/16+sin(x(:,2)/10);
                    l_Obj(:,2) = (x(:,1).^2+(x(:,2)-6).^4-2*x(:,1).*y(:,1)-(5-y(:,1)).^2)/80;          
            end
        end
        %% Calculate constraint violations
        function [u_Con,l_Con] = CalCon(obj,uDecs,lDecs,paratmeters,type)
            y = uDecs;
            x = lDecs;
            
            switch type
                case 'upper'
                    u_Con = (x(:,1)-0.5).^2+(x(:,2)-5).^2+(y(:,1)-5).^2-16;
                    l_Con = [];
                case 'lower'
                    u_Con = [];
                    l_Con(:,1) = -x(:,2)+x(:,1).^2;
                    l_Con(:,2) = 5*x(:,1).^2+x(:,2)-10;
                    l_Con(:,3) = x(:,2)+y(:,1)/6-5;
                    l_Con(:,4) = -x(:,1);
                case 'bilevel'
                    u_Con = -(1+x(:,1)+x(:,2));
                    l_Con = x(:,1).^2+ x(:,2).^2-y.^2;
            end
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,n)
            P = [];
        end
        
        function P = lower_PF(obj,y)
            P = [];
        end
    end
end