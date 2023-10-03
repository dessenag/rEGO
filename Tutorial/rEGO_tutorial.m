%% A tutorial on the refined Efficient Global Optimisation algorithm
% Please cite [1-3] when using this software for your work or research, Thank 
% you. 
% 
% This tutorial is the MATLAB code for "A global-local meta-modelling technique 
% for model updating" article [1].
%% The Test Functions
% The first test function is a modified version of the Branin function [4] with 
% an output scaled between 0 and 1 and two input variables. This function has 
% two local minima and its contour plot is shown in Fig. 1.
% 
% 
% 
% The second test function is the three storey structure from the EI at LANL 
% [5], the same as in [1]. Fig. 2 shows the system schematic and Fig. 3 the equivalent 
% mass-spring-damper system.
% 
% 
% 
% Fig. 2 Three storey structure schamatic and photos (Adapted from [5]).
% 
% 
% 
% Fig. 3 Equivalent Mass spring damper system.
% 
% For the sake of this tutorial, a reduction in stiffness of 25% in the second 
% interstory (${\mathit{\mathbf{k}}}_3$) is considered. Please note, only the 
% three stiffness values of the interstories (${\mathit{\mathbf{k}}}_{2-4}$) are 
% updated in this tutorial. The difference in the modal parameters between the 
% two systems, quantified as the modified total modal assurance criterion (MTMAC), 
% is to be minimised to characterise the change in stiffness. The numerical model 
% developed in [1] is used for this task.
%% The Modified Branin Function
% In this section the Modified Branin Function is minimised using the rEGO.

%% Preamble
close all
clear all
dd = split(fileparts(matlab.desktop.editor.getActiveFilename),'\Tutorial');
cd(dd{1}) %set working directory to main folder
addpath(genpath(dd{1})) % add paths to matlab only for this session
clear dd

%%-------------------------------------------------------------------------------------

fun = @(x)braninmodif_n(x); % assign Branin function
num_vari = 2; % call number of variables
lwb = [0 0]; upb = [1 1]; % search bounds
eps1 = 10^(-3); eps2 = 10^(-4); % stopping criterion 

[x,fval,results] = rEGO(fun,num_vari,lwb,upb,eps1,eps2)
%% The Three Storey Structure
% In this section the rEGO is used to identify the damage in a simulated structure 
% by minimising the MTMAC between the damaged and update system. Damage is simulated 
% with a 25% decrease in ${\mathit{\mathbf{k}}}_3$.
% 
% Let us load the dataset for the damaged mode and show the modal parameters 
% of the damaged and baseline model:

clear all
load LANL_3SS_dam_25_3.mat % load damaged modal parameters
load LANL_3SS.mat % load baseline modal parameters
%% 
% Modal properties of the baseline system:

baseline(:,2:end) 
%% 
% Modal properties of the damaged system:

damaged(:,2:end)
%% 
% Difference, in percentage, between damaged and baseline natural frequencies:

delta_w = 100.*(damaged(1,2:end)-baseline(1,2:end))./baseline(1,2:end)
%% 
% Diagonal of the MAC matrix of the damaged and baseline mode shapes:

mac = diag(compute_mac(damaged(2:end,2:end),baseline(2:end,2:end)))'
%% 
% $$\text{MTMAC}_{residuals}=1-\prod^n_{i=1}\frac{\text{MAC}(\bf{\phi}_i^E,\bf{\phi}_i^N)}{\Bigg(1+\frac{|\omega_i^N-\omega_i^E|}{|\omega_i^N+\omega_i^E|}\Bigg)}$$

yt = mac.*((1+abs((damaged(1,2:end)-baseline(1,2:end))./(damaged(1,2:end)+baseline(1,2:end)))).^(-1));
mtmac=1-prod(yt)
%% 
% By minimising the MTMAC, the damaged system can be identified starting from 
% the baseline system:

fun = @(x)frame_u_opti([ones(1,4) x],damaged(:,2:end),'mtmac'); % assign Branin function
num_vari = 3; % call number of variables
lwb = [.5 .5 .5]; upb = [1.01 1.01 1.01]; % search bounds
eps1 = 10^(-3); eps2 = 10^(-4); % stopping criterion 

[x,fval,results] = rEGO(fun,num_vari,lwb,upb,eps1,eps2)
% References
% [1] G. Dessena, D. I. Ignatyev, J. F. Whidborne, L. Zanotti Fragonara, A global-local 
% meta-modelling technique for model updating, Computer Methods in Applied Mechanics 
% and Engineering, Vol. (2023). (DOI: -)
% 
% [2] G. Dessena, D. I. Ignatyev, J. F. Whidborne, L. Zanotti Fragonara, A Kriging 
% Approach to Model Updating for Damage Detection, in European Workshop of Structural 
% Health Monitoring 2022, LNCE volume 254, EWSHM 2022, pp. 245–255 (2023). (DOI: 
% <https://doi.org/10.1007/978-3-031-07258-1_26 10.1007/978-3-031-07258-1_26>)
% 
% [3] G. Dessena, rEGO - A tutorial on the refined Efficient Global Optimisation, 
% GitHub (2023) (DOI: -)
% 
% [4] A. Forrester, A. Sobester, A. Keane, A., Engineering design via surrogate 
% modelling: a practical guide. Wiley. (DOI: <https://doi.org/10.1002/9780470770801 
% 10.1002/9780470770801>)
% 
% [5] E. Figueiredo, G. Park, J. Figueiras, C. Farrar, K. Worden, Structural 
% health monitoring algorithm comparisons using standard data sets, Los Alamos 
% National Laboratory (LANL), LA-14393 (DOI: <https://doi.org/10.2172/961604 10.2172/961604>)