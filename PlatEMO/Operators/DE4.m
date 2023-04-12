function Offspring = DE4(Parent1,Parent2,Parent3,Parent4,Parent5,type,Parameter)
%DE - Differential evolution operator.
%
%   Off = DE(P1,P2,P3) returns the offsprings generated by differential
%   evolution operator and polynomial mutation, where P1, P2, and P3 are
%   three sets of parents. If P1, P2, and P3 are arrays of INDIVIDUAL
%   objects, then Off is also an array of INDIVIDUAL objects; while if P1,
%   P2, and P3 are matrices of decision variables, then Off is also a
%   matrix of decision variables, i.e., the offsprings will not be
%   evaluated. Each object/row of P1, P2, and P3 is used to generate one
%   offspring by P1 + 0.5*(P2-P3).
%
%	Off = DE(P1,P2,P3,{CR,F,proM,disM}) specifies the parameters of
%	operators, where CR and F are the parameters in differental evolution,
%	proM is the expectation of number of bits doing mutation, and disM is
%	the distribution index of polynomial mutation.
%
%   Example:
%       Off = DE(Parent1,Parent2,Parent3)
%       Off = DE(Parent1.decs,Parent2.decs,Parent3.decs,{1,0.5,1,20})

%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 6
        [CR,F,proM,disM] = deal(Parameter{:});
    else       
        [CR,F,proM,disM] = deal(1,0.5,1,20);  
    end
    if isa(Parent1(1),'INDIVIDUAL')
        calObj  = true;
        Parent1 = Parent1.decs;
        Parent2 = Parent2.decs;
        Parent3 = Parent3.decs;
    else
        calObj = false;
    end
    [N,D]  = size(Parent1);
    Global = GLOBAL.GetObj();

    %% Differental evolution
    Site = rand(N,D) < CR;
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F/2*(Parent2(Site)-Parent3(Site))+ F/2*(Parent4(Site)-Parent5(Site));

    %% Polynomial mutation
    switch type
        case 'upper'
            lower = Global.u_lower;
            upper = Global.u_upper;
        case 'lower'
            lower = Global.l_lower;
            upper = Global.l_upper;
    end
    Lower = repmat(lower,N,1);
    Upper = repmat(upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    if calObj
        Offspring = INDIVIDUAL(Offspring);
    end
end