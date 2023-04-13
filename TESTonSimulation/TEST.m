function TEST(Global)
% <algorithm> <A>
% Approximation-guided evolutionary multi-objective algorithm II
% epsilon --- 0.1 --- The parameter in grid location calculation

%------------------------------- Reference --------------------------------
% M. Wagner and F. Neumann, A fast approximation-guided evolutionary
% multi-objective algorithm, Proceedings of the 15th Annual Conference on
% Genetic and Evolutionary Computation, 2013, 687-694.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    N = 1;
    X = round(rand(N,Global.u_D));
    Y = rand(N,Global.u_D,Global.l_D/Global.u_D);
    Y = round(reshape(Y.*X*40,N,Global.l_D));
    
      tic
    Pop= INDIVIDUAL(X,Y,'upper');
    toc
    tic
    Pop= INDIVIDUAL(Pop,'lower');
      toc
    
    tic
    Pop2 =INDIVIDUAL(X,Y,'bilevel');
    toc
    t =0;
end