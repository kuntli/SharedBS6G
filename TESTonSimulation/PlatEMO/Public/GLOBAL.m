classdef GLOBAL < handle
%GLOBAL - The class of experimental setting.
%
%   This is the class of experimental setting. An object of GLOBAL class
%   stores the algorithm to be executed, the problem to be solved, and all
%   the parameter settings like population size, number of objectives, and
%   maximum number of evaluations. This class also has several methods
%   which can be invoked by algorithms.
%
% GLOBAL properties:
%   Nu              <public>	upper level population size
%   u_M;                            % Number of upper level objectives
%   l_M;                            % Number of lower level objectives
%   u_D;                            % Number of upper level decision variables
%   l_D;                            % Number of lower level decision variables
%   u_lower         <read-only>	lower bound of each upper level decision variable
%   u_upper         <read-only>	upper bound of each upper level decision variable
%   l_lower         <read-only>	lower bound of each lower level decision variable
%   l_upper         <read-only>	upper bound of each lower level decision variable
%   algorithm       <read-only>	algorithm function
%   problem         <read-only>	problem function
%   encoding        <read-only> encoding of the problem
%   u_evaluated     <read-only>	number of evaluated upper level individuals
%   l_evaluated     <read-only>	number of evaluated lower level individuals
%   u_evaluation    <read-only>	maximum number of upper level evaluations
%   l_evaluation    <read-only>	maximum number of lower level evaluations
%   gen             <read-only>	current generation
%   maxgen          <read-only>	maximum generation
%   run             <read-only>	run number
%   runtime         <read-only>	runtime
%   save            <read-only> number of saved populations
%   result          <read-only>	set of saved populations
%   PF              <read-only>	true Pareto front
%   parameter       <read-only>	parameters of functions specified by users
%   outputFcn    	<read-only>	function invoked after each generation
%
% GLOBAL methods:
%   GLOBAL          <public>    the constructor, all the properties will be
%                               set when the object is creating
%   Start           <public>    start running the algorithm
%   Initialization  <public>    randomly generate an initial population
%   NotTermination  <public>	terminate the algorithm if the number of
%                               evaluations has exceeded
%   ParameterSet    <public>    obtain the parameter settings from user
%   GetObj          <static>    get the current GLOBAL object

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        Nu          = 50;               % upper level Population size
        Nl          =  50;               % lower level Population size
    end
    properties(SetAccess = ?PROBLEM)
        u_M;                            % Number of upper level objectives
        l_M;                            % Number of lower level objectives
        u_D;                            % Number of upper level decision variables
        l_D;                            % Number of lower level decision variables
        u_lower;                        % Lower bound of each upper level decision variable
        u_upper;                        % Upper bound of each upper level decision variable
        l_lower;                        % Lower bound of each lower level decision variable
        l_upper;                        % Upper bound of each lower level decision variable
        encoding   = 'real';            % Encoding of the problem
        u_evaluation = 1000000000;           % Maximum number of upper level evaluations
        l_evaluation = 1000000000;        % Maximum number of lower level evaluations
    end
    properties(SetAccess = ?INDIVIDUAL)
        u_evaluated  = 0;               % Number of upper level evaluated individuals
        l_evaluated  = 0;               % Number of lower level evaluated individuals
    end
    properties(SetAccess = private)
        algorithm  = @NSGAII;       	% Algorithm function
        problem    = @DTLZ2;            % Problem function
        gen;                            % Current generation
        maxgen;                         % Maximum generation
        run        = 1;                 % Run number
        runtime    = 0;                 % Runtime
        save       = 0;             	% Number of saved populations
        result     = {};                % Set of saved populations
        PF;                             % True Pareto front
        parameter  = struct();      	% Parameters of functions specified by users
        outputFcn  = @GLOBAL.Output;  	% Function invoked after each generation
    end
    methods
        %% Constructor
        function obj = GLOBAL(varargin)
        %GLOBAL - Constructor of GLOBAL class.
        %
        %   GLOBAL('-Name',Value,'-Name',Value,...) returns an object with
        %   the properties specified by the inputs.
        %
        %   Example:
        %       GLOBAL('-algorithm',@NSGAII,'-problem',@DTLZ2,'-N',100,...
        %              '-M',2,'-D',10)
        
            obj.GetObj(obj);
            % Initialize the parameters which can be specified by users
            propertyStr = {'Nu','Nl','u_M','l_M','u_D','l_D','algorithm','problem','u_evaluation','l_evaluation','run','save','outputFcn'};
            if nargin > 0
                IsString = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
                [~,Loc]  = ismember(varargin(IsString),cellfun(@(S)['-',S],propertyStr,'UniformOutput',false));
                for i = find(Loc)
                    obj.(varargin{IsString(i)}(2:end)) = varargin{IsString(i)+1};
                end
            end
            % Instantiate a problem object
            obj.problem = obj.problem();
            % Add the folders of the algorithm and problem to the top of
            % the search path
            addpath(fileparts(which(class(obj.problem))));
            addpath(fileparts(which(func2str(obj.algorithm))));
        end
        %% Start running the algorithm
        function Start(obj)
        %Start - Start running the algorithm.
        %
        %   obj.Start() runs the algorithm. This method of one GLOBAL
        %   object can only be invoked once.
        %
        %   Example:
        %       obj.Start()

            if obj.u_evaluated <= 0
                obj.PF = obj.problem.PF();
                try
                    tic;
                    obj.algorithm(obj);
                catch err
                    if strcmp(err.identifier,'GLOBAL:Termination')
                        return;
                    else
                        rethrow(err);
                    end
                end
                obj.u_evaluated = max(obj.u_evaluated,obj.u_evaluation);
                if isempty(obj.result)
                    obj.result = {obj.u_evaluated,obj.l_evaluated,INDIVIDUAL()};
                end
            	obj.outputFcn(obj);
            end
        end
        %% Randomly generate an initial population
        function Population = Initialization(obj,N,type)
        %Initialization - Randomly generate an initial population.
        %
        %   P = obj.Initialization() returns an initial population, i.e.,
        %   an array of obj.N INDIVIDUAL objects.
        %
        %   P = obj.Initialization(N) returns an initial population with N
        %   individuals.
        %
        %   Example:
        %       P = obj.Initialization()
        
            if nargin < 2
                N = obj.N;
            end
            Population = INDIVIDUAL(obj.problem.Init(N,type));
        end
        %% Terminate the algorithm if the number of evaluations has exceeded
        function notermination = NotTermination(obj,Output,notermination)
        %NotTermination - Terminate the algorithm if the number of
        %evaluations has exceeded.
        %
        %   obj.NotTermination(P) stores the population P as the final
        %   population for output, and returns true if the number of
        %   evaluations has not exceeded (otherwise returns false). If the
        %   number of evaluations has exceeded, then throw an error to
        %   terminate the algorithm forcibly.
        %
        %   This function should be invoked after each generation, and the
        %   function obj.outputFcn will be invoked when obj.NotTermination
        %   is invoked.
        %
        %   The runtime of executing this function and obj.outputFcn is not
        %   counted as part of the runtime of the algorithm.
        %
        %   Example:
        %       obj.NotTermination(Population)
        
            % Accumulate the runtime
            obj.runtime = obj.runtime + toc;
            % Save the last population
            if obj.save<=0; num=10; else; num=obj.save; end
            index = max(1,min(min(num,size(obj.result,1)+1),ceil(num*obj.u_evaluated/obj.u_evaluation)));
            obj.result(index,:) = {obj.u_evaluated,obj.l_evaluated,Output,obj.runtime};
            % Invoke obj.outputFcn
            drawnow();
            obj.outputFcn(obj);
            % Detect whether the number of evaluations has exceeded
%            notermination = obj.evaluated < obj.evaluation;
%            assert(notermination,'GLOBAL:Termination','Algorithm has terminated');
            tic;
        end
        %% Obtain the parameter settings from user
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameter settings from user.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) returns the values of
        %   p1, p2, ..., where v1, v2, ... are their default values. The
        %   values are specified by the user with the following form:
        %   MOEA(...,{@X,p1,p2,...},...), where X is the function name of
        %   the caller.
        %
        %   MOEA(...,{@X,[],p2,...},...) indicates that p1 is not specified
        %   by the user, and p1 equals to its default value v1.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            CallStack = dbstack();
            caller    = CallStack(2).file;
            caller    = caller(1:end-2);
            varargout = varargin;
            if isfield(obj.parameter,caller)
                specified = ~cellfun(@isempty,obj.parameter.(caller));
                varargout(specified) = obj.parameter.(caller)(specified);
            end
        end
        %% Variable constraint
        function set.Nu(obj,value)
            obj.Validation(value,'int','size of upper level population ''-Nu''',1);
            obj.Nu = value;
        end
        function set.Nl(obj,value)
            obj.Validation(value,'int','size of lower level population ''-Nl''',1);
            obj.Nl = value;
        end
        function set.u_M(obj,value)
            obj.Validation(value,'int','number of upper level objectives ''-u_M''',2);
            obj.u_M = value;
        end
        function set.l_M(obj,value)
            obj.Validation(value,'int','number of lower level objectives ''-l_M''',2);
            obj.l_M = value;
        end
        function set.u_D(obj,value)
            obj.Validation(value,'int','number of upper level variables ''-u_D''',1);
            obj.u_D = value;
        end
        function set.l_D(obj,value)
            obj.Validation(value,'int','number of lower level variables ''-l_D''',1);
            obj.l_D = value;
        end
        function set.algorithm(obj,value)
            if iscell(value) 
                obj.Validation(value{1},'function','algorithm ''-algorithm''');
                obj.algorithm = value{1};
                obj.parameter.(func2str(value{1})) = value(2:end);
            else
                obj.Validation(value,'function','algorithm ''-algorithm''');
                obj.algorithm = value;
            end
        end
        function set.problem(obj,value)
            if iscell(value)
                obj.Validation(value{1},'function','test problem ''-problem''');
                obj.problem = value{1};
                obj.parameter.(func2str(value{1})) = value(2:end);
            elseif ~isa(value,'PROBLEM')
                obj.Validation(value,'function','test problem ''-problem''');
                obj.problem = value;
            else
                obj.problem = value;
            end
        end
        function set.u_evaluation(obj,value)
            obj.Validation(value,'int','number of upper level evaluations ''u_evaluation''',1);
            obj.u_evaluation = value;
        end
        function set.l_evaluation(obj,value)
            obj.Validation(value,'int','number of lower level evaluations ''-l_evaluation''',1);
            obj.l_evaluation = value;
        end
        function set.run(obj,value)
            obj.Validation(value,'int','run number ''-run''',1);
            obj.run = value;
        end
        function set.save(obj,value)
            obj.Validation(value,'int','number of saved populations ''-save''',0);
            obj.save = value;
        end
        %% Variable dependence
        function value = get.gen(obj)
            value = ceil(obj.u_evaluated/obj.N);
        end
        function value = get.maxgen(obj)
            value = ceil(obj.u_evaluation/obj.N);
        end
    end
    methods(Static)
        %% Get the current GLOBAL object
        function obj = GetObj(obj)
        %GetObj - Get the current GLOBAL object.
        %
        %   Global = GLOBAL.GetObj() returns the current GLOBAL object.
        %
        %   Example:
        %       Global = GLOBAL.GetObj()
        
            persistent Global;
            if nargin > 0
                Global = obj;
            else
                obj = Global;
            end
        end
    end

    % The following functions cannot be invoked by users
    methods(Access = private)
        %% Check the validity of the specific variable
        function Validation(obj,value,Type,str,varargin)
            switch Type
                case 'function'
                    assert(isa(value,'function_handle'),'INPUT ERROR: the %s must be a function handle',str);
                    assert(~isempty(which(func2str(value))),'INPUT ERROR: the function <%s> does not exist',func2str(value));
                case 'int'
                    assert(isa(value,'double') && isreal(value) && isscalar(value) && value==fix(value),'INPUT ERROR: the %s must be an integer scalar',str);
                    if ~isempty(varargin); assert(value>=varargin{1},'INPUT ERROR: the %s must be not less than %d',str,varargin{1}); end
                    if length(varargin) > 1; assert(value<=varargin{2},'INPUT ERROR: the %s must be not more than %d',str,varargin{2}); end
                    if length(varargin) > 2; assert(mod(value,varargin{3})==0,'INPUT ERROR: the %s must be a multiple of %d',str,varargin{3}); end
            end
        end
    end
    methods(Access = private, Static)
        %% Display or save the result after the algorithm is terminated
        function Output(obj)
            clc; fprintf('%s on %s, %d upper level objectives %d lower level objectives, %d upper level variables %d lower level variables, run %d , %.2fs passed...\n',...
                         func2str(obj.algorithm),class(obj.problem),obj.u_M,obj.l_M,obj.u_D,obj.l_D,obj.run,obj.runtime);
            if obj.u_evaluated >= obj.u_evaluation
                if obj.save == 0
                    % Identify the feasible and non-dominated solutions in the
                    % final population
                    Feasible     = find(all(obj.result{end}.ucons<=0,2));
                    NonDominated = NDSort(obj.result{end}(Feasible).uobjs,1) == 1;
                    Population   = obj.result{end}(Feasible(NonDominated));
                    % Calculate the metric values
                    if length(Population) >= size(obj.PF,1)
                        Metrics = {@HV};
                    else
                        Metrics = {@IGD};
                    end
                    Score     = cellfun(@(S)GLOBAL.Metric(S,Population,obj.PF),Metrics,'UniformOutput',false);
                    MetricStr = cellfun(@(S)[func2str(S),' : %.4e  '],Metrics,'UniformOutput',false);
                    % Display the results
                    figure('NumberTitle','off','UserData',struct(),...
                           'Name',sprintf([MetricStr{:},'Runtime : %.2fs'],Score{:},obj.runtime));
                    title(sprintf('%s on %s',func2str(obj.algorithm),class(obj.problem)),'Interpreter','none');
                    Draw(Population.uobjs);
                    % Add new menus to the figure
                    top = uimenu(gcf,'Label','Data Source');
                    uimenu(top,'Label','Result (PF)',     'CallBack',{@(hObject,~,obj,P)eval('cla;Draw(P.objs);GLOBAL.cb_menu(hObject);'),obj,Population},'Checked','on');
                    uimenu(top,'Label','Result (PS)',     'CallBack',{@(hObject,~,obj,P)eval('cla;Draw(P.decs);GLOBAL.cb_menu(hObject);'),obj,Population});
                    uimenu(top,'Label','Result (Special)','CallBack',{@(hObject,~,obj,P)eval('obj.problem.Draw(P.decs);GLOBAL.cb_menu(hObject);'),obj,Population});
                    uimenu(top,'Label','True PF',         'CallBack',{@(hObject,~,obj)eval('cla;Draw(obj.PF);GLOBAL.cb_menu(hObject);'),obj},'Separator','on');
                    uimenu(top,'Label','IGD',             'CallBack',{@GLOBAL.cb_metric,obj,@IGD},'Separator','on');
                    uimenu(top,'Label','HV',              'CallBack',{@GLOBAL.cb_metric,obj,@HV});
                    uimenu(top,'Label','GD',              'CallBack',{@GLOBAL.cb_metric,obj,@GD});
                    uimenu(top,'Label','CPF',             'CallBack',{@GLOBAL.cb_metric,obj,@CPF});
                else
                    folder = fullfile('Data',func2str(obj.algorithm));
                    [~,~]  = mkdir(folder);
                    result         = obj.result;
                    metric.runtime = obj.runtime;
                    save(fullfile(folder,sprintf('%s_%s_uM%d_lM%d_uD%d_lD%d_%d.mat',func2str(obj.algorithm),class(obj.problem),obj.u_M,obj.l_M,obj.u_D,obj.l_D,obj.run)),'result','metric');
                end
            end
        end
        function cb_metric(hObject,eventdata,obj,metric)
            metricName   = func2str(metric);
            MetricValues = get(gcbf,'UserData');
            % Calculate the specified metric value of each population
            if ~isfield(MetricValues,func2str(metric)) 
                tempText = text('Units','normalized','Position',[.4 .5 0],'String','Please wait ... ...'); drawnow();
                MetricValues.(metricName)(:,1) = obj.result(:,1);
                MetricValues.(metricName)(:,2) = cellfun(@(S)GLOBAL.Metric(metric,S,obj.PF),obj.result(:,2),'UniformOutput',false);
                set(gcbf,'UserData',MetricValues);
                delete(tempText);
            end
            % Display the results
            cla; Draw(cell2mat(MetricValues.(metricName)),'-k.','LineWidth',1.5,'MarkerSize',10);
            xlabel('Number of Evaluations');
            ylabel(metricName);
            GLOBAL.cb_menu(hObject);
        end
        function cb_menu(hObject)
            % Switch the selected menu
            set(get(get(hObject,'Parent'),'Children'),'Checked','off');
            set(hObject,'Checked','on');
        end
        function value = Metric(metric,Population,PF)
            % Calculate the metric value of the population
            Feasible     = find(all(Population.cons<=0,2));
            NonDominated = NDSort(Population(Feasible).objs,1) == 1;
            try
                value = metric(Population(Feasible(NonDominated)).objs,PF);
            catch
                value = NaN;
            end
        end
    end
end