classdef INDIVIDUAL < handle
%INDIVIDUAL - The class of an individual.
%
%   This is the class of an individual. An object of INDIVIDUAL stores all
%   the properties including decision variables, objective values,
%   constraint violations, and additional properties of an individual.
%
% INDIVIDUAL properties:
%   udec         <read-only>    upper level decision variables of the individual
%   ldec         <read-only>    lower level decision variables of the individual
%   uobj         <read-only>    upper level objective values of the individual
%   lobj         <read-only>    lower level objective values of the individual
%   ucon         <read-only>    upper level constraint violations of the individual
%   lcon         <read-only>    lower level constraint violations of the individual
%   uFNo         <read-only>    upper level nondominated sorting rank of the individual
%   lFNo         <read-only>    lower level nondominated sorting rank of the individual
%   add         <read-only>     additional properties of the individual
%
% INDIVIDUAL methods:
%   INDIVIDUAL	<public>        the constructor, all the properties will be
%                               set when the object is creating
%   udecs        <public>      	get the matrix of upper level decision variables of the
%                               population
%   ldecs        <public>      	get the matrix of lower level decision variables of the
%                               population
%   uobjs        <public>       get the matrix of upper level objective values of the
%                               population
%   lobjs        <public>       get the matrix of lower level objective values of the
%                               population
%   ucons        <public>       get the matrix of upper level constraint violations of
%                               the population
%   lcons        <public>       get the matrix of lower level constraint violations of
%                               the population
%   uFNos        <public>       get the vector of upper level nondominated sorting rank of
%                               the population
%   lFNos        <public>       get the vector of lower level nondominated sorting rank of
%                               the population
%   adds        <public>        get the matrix of additional properties of
%                               the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        udec;       % upper level Decision variables of the individual
        ldec;       % lower level Decision variables of the individual
        uobj;       % upper level Objective values of the individual
        lobj;       % lower level Objective values of the individual
        ucon;       % upper level Constraint violations of the individual
        lcon;       % lower level Constraint violations of the individual
        uFNo        % upper level nondominated sorting rank of the individual
        lFNo        % lower level nondominated sorting rank of the individual
        parameter;
        add;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = INDIVIDUAL(varargin)
            %INDIVIDUAL - Constructor of INDIVIDUAL class.
            %
            %   H = INDIVIDUAL(Dec) creates an array of individuals (i.e., a
            %   population), where Dec is the matrix of decision variables of
            %   the population. The objective values and constraint violations
            %   are automatically calculated by the test problem functions.
            %   After creating the individuals, the number of evaluations will
            %   be increased by length(H).
            %
            %   H = INDIVIDUAL(Dec,AddProper) creates the population with
            %   additional properties stored in AddProper, such as the velocity
            %   in particle swarm optimization.
            %
            %   Example:
            %       H = INDIVIDUAL(rand(100,3))
            %       H = INDIVIDUAL(rand(100,10),randn(100,3))
            if nargin > 0
                if ~isa(varargin{1}(1),'INDIVIDUAL')
                    uDecs = varargin{1};
                    lDecs = varargin{2};
                    type  = varargin{3};
                    
                    % Create new objects
                    obj(1,size(uDecs,1)) = INDIVIDUAL;
                    Global = GLOBAL.GetObj();
                    % Set the infeasible decision variables to boundary values
                    if ~isempty(Global.u_lower) && ~isempty(Global.u_upper)
                        u_Lower = repmat(Global.u_lower,length(obj),1);
                        u_Upper = repmat(Global.u_upper,length(obj),1);
                        uDecs  = max(min(uDecs,u_Upper),u_Lower);
                    end
                    if ~isempty(Global.l_lower) && ~isempty(Global.l_upper)
                        l_Lower = repmat(Global.l_lower,length(obj),1);
                        l_Upper = repmat(Global.l_upper,length(obj),1);
                        lDecs  = max(min(lDecs,l_Upper),l_Lower);
                    end
                    % Calculte the objective values and constraint violations
                    %                Decs = Global.problem.CalDec(Decs);
                    [u_Objs,l_Objs,parameters]= Global.problem.CalObj(uDecs,lDecs,obj.parameters,type);
                    [u_Cons,l_Cons] = Global.problem.CalCon(uDecs,lDecs,parameters,type);
                    % Assign the decision variables, objective values,
                    % constraint violations, and additional properties
                    for i = 1 : length(obj)
                        obj(i).udec = uDecs(i,:);
                        obj(i).ldec = lDecs(i,:);
                        switch type
                            case 'upper'
                                obj(i).uobj = u_Objs(i,:);
                                obj(i).ucon = u_Cons(i,:);
                                if ~isempty(parameters)
                                obj(i).parameter = parameters(i);
                                end
                            case 'lower'
                                obj(i).lobj = l_Objs(i,:);
                                obj(i).lcon = l_Cons(i,:);
                                if ~isempty(parameters)
                                obj(i).parameter = parameters(i);
                                end
                            case 'bilevel'
                                obj(i).uobj = u_Objs(i,:);
                                obj(i).ucon = u_Cons(i,:);
                                
                                obj(i).lobj = l_Objs(i,:);
                                obj(i).lcon = l_Cons(i,:);
                                if ~isempty(parameters)
                                obj(i).parameter = parameters(i);
                                end
                        end
                    end
                    if nargin > 3
                        AddProper = varargin{4};
                        for i = 1 : length(obj)
                            obj(i).add = AddProper(i,:);
                        end
                    end
                    % Increase the number of evaluated individuals
                    switch type
                        case 'upper'
                            Global.u_evaluated = Global.u_evaluated + length(obj);
                        case 'lower'
                            Global.l_evaluated = Global.l_evaluated + length(obj);
                        case 'bilevel'
                            Global.u_evaluated = Global.u_evaluated + length(obj);
                            Global.l_evaluated = Global.l_evaluated + length(obj);
                    end
                else
                    obj = varargin{1};
                    type = varargin{2};
                    Global = GLOBAL.GetObj();
                    uDecs = obj.udecs;
                    lDecs = obj.ldecs;
                    parameters = obj.parameters;
                    
                    [u_Objs,l_Objs,parameters] = Global.problem.CalObj(uDecs,lDecs,parameters,type);
                    [u_Cons,l_Cons] = Global.problem.CalCon(uDecs,lDecs,parameters,type);
                    for i = 1 : length(obj)
                        switch type
                            case 'upper'
                                obj(i).uobj = u_Objs(i,:);
                                obj(i).ucon = u_Cons(i,:);
                                if ~isempty(parameters)
                                obj(i).parameter = parameters(i);   
                                end
                            case 'lower'
                                obj(i).lobj = l_Objs(i,:);
                                obj(i).lcon = l_Cons(i,:);
                                if ~isempty(parameters)
                                obj(i).parameter = parameters(i);  
                                end
                        end
                    end
                    % Increase the number of evaluated individuals
                    switch type
                        case 'upper'
                            Global.u_evaluated = Global.u_evaluated + length(obj);
                        case 'lower'
                            Global.l_evaluated = Global.l_evaluated + length(obj);
                    end
                end
            end
        end
        
        
        function obj = add_FNo(obj,FNo,type)          
            % Increase the number of evaluated individuals
            switch type
                case 'upper'
                    for i = 1:length(obj)
                       obj(i).uFNo = FNo(i); 
                    end
                case 'lower'
                    for i = 1:length(obj)
                       obj(i).lFNo = FNo(i); 
                    end
            end
        end
        %% Get the matrix of upper level decision variables of the population
        function value = udecs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.udec);
        end
         %% Get the matrix of lower level decision variables of the population
        function value = ldecs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.ldec);
        end
        %% Get the matrix of upper level objective values of the population
        function value = uobjs(obj)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.uobj);
        end
          %% Get the matrix of lower level objective values of the population
        function value = lobjs(obj)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.lobj);
        end
        %% Get the matrix of upper level constraint violations of the population
        function value = ucons(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.ucon);
        end
          %% Get the matrix of lower level constraint violations of the population
        function value = lcons(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.lcon);
        end
          %% Get the vector of upper level nondominated sorting rank of the population
        function value = uFNos(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.uFNo);
        end
          %% Get the vector of lower level nondominated sorting rank of the population
        function value = lFNos(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.lFNo);
        end    
        %%
        function value = parameters(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(2,obj.parameter);
        end 
        %% Get the matrix of additional properties of the population
        function [value,obj] = adds(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.
            if nargin > 1
                for i = 1 : length(obj)               
                    obj(i).add = AddProper(i,:);
                end
            end
            value = cat(1,obj.add);
        end
    end
end