function [x,fval,results] = rEGO(fun_name,num_vari,lwb,upb,eps1,eps2,printing)
%% [x,fval,results] = rEGO(fun_name,num_vari,lwb,upb,eps,printing)

%% Please cite the works under "References" when using this program.

% The program computes the modal properties, the state space matrices, and 
% the fitted model of the given FRFs via the Loewner Framework. The state
% space approximation has the form: H(s)=C(sE-A)^(-1)B.
%     
%     Given: 
%     fun_name = function of which we are minimising - must be a single-
%                objective problem;
%     num_vari = number of variables of the problem;
%     lwb      = lower bound of the problem;
%     upb      = upper bound of the problem;
%     eps1     = first stopping criterion (optional - 0.01);
%     eps2     = second stopping criterion (optional - 10^(-4));
%     printing = (optional - set as 1 as standard)if = 1 prints and displys
%     information during the terations;
%
%     Returns:
%     x:       position, in terms of variables, of the minimum;
%     fval:    value of the minimum;
%     results: initial: set of points from the design of experiment
%              samples: values of the function from DoE
%              min_x: values of the variables for the minimum objective at
%              each iteration
%              min_y: minimum value of the function at each iteration
%              x: last set of variables in the data pool
%              y: last set of function values in the data pool
%              evaluation: number of function evaluations to convergence
%              fval: value of the minimum;
%              iteration: number of iterations to convergence
%              EI: expected improvement at convergence
%              diff: cartesian distance between subsequent minima at
%              convergence
%              flag: if 1, iteration stopped because 100*number of variables
%              without improvement
%              refinement: number of search space refinements
%              stall: number of iteration without an improvement

%% Disclaimer
% This program is free software: you can redistribute it and/or modify  it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License and GNU
% Lesser General Public License along with this program. If not, see
% <http://www.gnu.org/licenses/>.

%% Credits
% Implementation of the refined Efficient Global Optimisation:
% Gabriele Dessena
% gdessena@ing.uc3m.es
% Universidad Carlos III de Madrid
% 2023/10/02 v0.1 - pre release
% 
% The code is based on:
% - The repository "Single_objective_EGO_algorithms"
%   algorithms by Qi Zhang
%   (https://github.com/109902249/Single_objective_EGO_algorithms)
% - The dace Toolbox by Hans Bruun Nielsen, Søren Nymand and 
%   Lophaven Jacob Søndergaard
%   (https://github.com/psbiomech/dace-toolbox-source)
%
%% Changelog
% Reserved
%% References
%% Please cite the works under "References" when using this program
%  [1] G. Dessena, D. I. Ignatyev, J. F. Whidborne, L. Zanotti Fragonara, A
%      global-local meta-modelling technique for model updating,
%      Computer Methods in Applied Mechanics and Engineering, Vol. 
%      (2023). (DOI: -)
%  [2] G. Dessena, D. I. Ignatyev, J. F. Whidborne, L. Zanotti Fragonara, 
%      A Kriging Approach to Model Updating for Damage Detection, LNCE 254,
%      EWSHM 2022, pp. 245-255 (2023). (DOI: 10.1007/978-3-031-07258-1_26)
%  [3] G. Dessena, rEGO - A tutorial for the refined Efficient Optimisation 
%      algorithm, Git Hub,(2023).
%      (DOI: )
%% Please cite the works under "References" when using this program
%--------------------------------------------------------------------------

if ~exist('eps1','var')
    % if parameter does not exist, so default it to something
    eps1 = 0.01;
end
if ~exist('eps2','var')
    % if parameter does not exist, so default it to something
    eps2 = 10e-4;
end
if ~exist('printing','var')
    % if parameter does not exist, so default it to something
    p_flag = 1;
end
%--------------------------------------------------------------------------
%set up design space
design_space=[lwb;upb];
dso = design_space;
optimum = 0; %stopping flag
%--------------------------------------------------------------------------
% the number of initial design points
num_initial = num_vari*10;
% variable for initialisaion of vector f_min
max_evaluation = 100;
%--------------------------------------------------------------------------
% initial design points using Latin hypercube sampling method
sample_x = repmat(design_space(1,:),num_initial,1) + ...
    repmat(design_space(2,:)-design_space(1,:),num_initial,1).*...
    lhsdesign(num_initial,num_vari,'criterion','maximin',...
    'iterations',1000);
[m,~] = size(sample_x);
sample_y = zeros(m,1); %create empty vector for y=f(x)
for i = 1:m
    sample_y(i,1) = feval(fun_name,sample_x(i,:));
end
results.initial = [sample_x sample_y];
% record the f_min in each iteration
f_min = zeros(max_evaluation - num_initial + 1,1);
% the current best solution
f_min(1) = min(sample_y);
% the current iteration and evaluation
evaluation = size(sample_x,1);
% init stopping flags
iteration = 0;
refinement = 0;
stall = 0;
stagn = 100*num_vari; % stagnation condition
%--------------------------------------------------------------------------
if p_flag == 1
    % plot the iteration history
    plot(iteration,f_min(1),'r-o');title(sprintf('iteration: %d, evaluations:%d',iteration,evaluation));drawnow;
    % print the current information to the screen
    fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(1), optimum);
end
%--------------------------------------------------------------------------
% record progresses
[~,ijk] = min(sample_y);
results.min_x(iteration+1,:) = sample_x(ijk,:);
results.min_y(iteration+1)= min(sample_y);
%--------------------------------------------------------------------------
% the iteration procedure
%initialise iteration flags
stopping = 1;
count_m = 0;
count_r = 0;
GAOptions = optimoptions('ga','FunctionTolerance',1e-9,...
    'PopulationSize',10*length(lwb),'Display','off',...
    'UseParallel',false); %init GA parameters

while stopping >  0 % verify convergence not reached
    if iteration > 0 % keeps track of iteration number
        if abs(my - min(sample_y))<10e-10 % verify stall iteration
            stall = stall + 1;
            if stall > stagn-1 % stagnation stopping condition
                % saving results
                results.x = sample_x;
                results.y = sample_y;
                results.evaluation = evaluation;
                [~,ijk] = min(sample_y);
                results.min_x(end+1,:) = sample_x(ijk,:);
                results.min_y(end+1)= min(sample_y);
                results.fval=fval;
                results.iteration = iteration;
                results.EI = EI;
                results.flag = 0;
                results.refinement = refinement;
                results.stall = 0;
                stopping = 0;
            end
        else
            stall = 1;
        end
        if exist('miny','var') == 1 %check for the conservation of the minimum
            if sum(sample_y==miny)<1
                sample_y=[sample_y;miny];
                sample_x=[sample_x;minx];
            end
        end
        [sample_x,ia] = unique(sample_x,'rows'); %check for uniqueness
        sample_y=sample_y(ia,:);
        % build (or rebuild) the initial Kriging model
        [kriging_model, ~] = dacefit(sample_x,sample_y,'regpoly0','corrgauss',ones(1,num_vari),(10^(-3))*ones(1,num_vari),(10^2)*ones(1,num_vari));

    else %first iteration condition
        my = min(sample_y);
        meno = my;
        % build (or rebuild) the initial Kriging model
        [kriging_model, ~] = dacefit(sample_x,sample_y,'regpoly0','corrgauss',ones(1,num_vari),(10^(-3))*ones(1,num_vari),(10^2)*ones(1,num_vari));
    end
    iteration = iteration + 1; %counting iterations
    results.samples(iteration) = length(sample_y);

    % the Expected Improvement criterion
    infill_criterion = @(x)Infill_Standard_EI(x,kriging_model,min(sample_y));

    % find the point with the highest EI value using GA algorithm
    [best_x, EI]=ga(infill_criterion,num_vari,[],[],[],[],...
        design_space(1,:),design_space(2,:),[],GAOptions);
    disp('EI: '+string(max(abs(EI))))

    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    evaluation = evaluation+1;

    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    f_min(iteration+1) = min(sample_y);
    if p_flag == 1
        % plot the iteration history
        plot(0:iteration,f_min(1:iteration+1),'r-o');title(sprintf('iteration: %d, evaluations:%d',iteration,evaluation));drawnow;
        % print the current information to the screen
        fprintf(' iteration: %d, evaluation: %d, current best solution: %f, real optimum: %f\n', iteration, evaluation, f_min(iteration+1), optimum);
    end

    if max(abs(EI)) < eps1 || count_r >= 100*num_vari %First Stopping Criterion
        if count_r >= 100*num_vari
            disp("Refinement by max iterations without refinements")
            count_r = 0;
        end
        % enforce uniqueness
        [sample_x,ia] = unique(sample_x,'rows');
        sample_y=sample_y(ia,:);
        % infill Kriging model
        [kriging_model, ~] = dacefit(sample_x,sample_y,'regpoly0',...
            'corrgauss',ones(1,num_vari),(10^(-3))*ones(1,num_vari),...
            (10^2)*ones(1,num_vari));
        % predict minimum
        f_pred = @(x)predictor(x, kriging_model);
        [pred_x, dummi]=ga(f_pred,num_vari,[],[],[],[],...
            design_space(1,:),design_space(2,:),[],GAOptions);
        % add points to set
        sample_x = [sample_x;pred_x];
        pred_y = feval(fun_name,pred_x);
        [~,idx_m] = min(sample_y);
        sample_y = [sample_y;pred_y];
        % update counting
        evaluation = evaluation+1;
        f_min(iteration+1) = min(sample_y);
        results.min_x(iteration+1)=f_min(end);
        % compute local search
        diff = pdist([sample_x(idx_m,:);pred_x],'euclidean');
        if p_flag == 1
            fprintf(' iteration: %d, evaluation: %d, current best solution: %f, difference: %f\n, samples: %f\n', iteration, evaluation, f_min(iteration+1), abs(pred_y-dummi),length(sample_y));
            disp('Refinement:'+string(refinement))
            disp('Distance is ' + string(diff))
        end
        if diff<=eps2% second stopping criterion
            results.x = sample_x;
            results.y = sample_y;
            results.evaluation = evaluation;
            [fval,id] = min(sample_y);
            x = sample_x(id,:);
            results.fval=fval;
            results.iteration = iteration;
            results.EI = EI;
            results.diff = diff;
            results.flag = 0;
            results.refinement = refinement;
            results.stall = 0;
            [results.min_y(end+1),ijk]=min(sample_y);
            results.min_x(end+1,:) = sample_x(ijk,:);
            stopping = 0;
        elseif length(sample_y) > num_initial-1 % search space refinement
            refinement = refinement+1;
            % detect minimum
            [miny,b]=min(sample_y);
            minx=sample_x(b,:);
            % refine search space
            ub = minx+.25.*(design_space(2,:)-design_space(1,:));
            lb = minx-.25.*(design_space(2,:)-design_space(1,:));
            ub(ub>dso(2,:)) = dso(2,ub>dso(2,:));
            lb(lb<dso(1,:)) = dso(1,lb<dso(1,:));
            design_space = [lb;ub];
            % find points in search space
            k = find(sum([sample_x < lb sample_x > ub],2)==0);
            sample_x = sample_x(k,:);
            sample_y = sample_y(k,:);
            % select points in differents blocks
            if length(k) >=75*num_vari-1
                x = ga(@(x)dsobj(x,sample_x, sample_y,25*num_vari),1,...
                    [],[],[],[],[10^(-50)],[.5],[],GAOptions);
                [sample_x, sample_y] = dsmerge(sample_x, sample_y, x, 2, 3);
                clear k
            end
            % ensure minimun values stays
            sample_y=[sample_y;miny];
            sample_x=[sample_x;minx];
        end
    end
    % record progresses
    [results.min_y(iteration+1),ijk] = min(sample_y);
    results.min_x(iteration+1,:) = sample_x(ijk,:);
    if meno == min(sample_y)
        count_m = 1 +  count_m;
        disp("Stall iterations "+string(count_m))
    else
        count_m = 1;
    end
    meno = min(sample_y);
    if count_m > 100*num_vari
        results.x = sample_x;
        results.y = sample_y;
        results.evaluation = evaluation;
        [fval,id] = min(sample_y);
        x = sample_x(id,:);
        results.fval=fval;
        results.iteration = iteration;
        results.EI = EI;
        %        results.diff = diff;
        results.flag = 1;
        results.refinement = refinement;
        results.stall = 1;

        disp('Maximum stall generation 100 x n_variables')
        stopping = 0;
    end
    count_r = count_r +1;
end

end


