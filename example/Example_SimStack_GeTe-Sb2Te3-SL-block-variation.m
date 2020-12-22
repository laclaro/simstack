 %{ 
 %  This file is part of the SimStack distribution (https://github.com/laclaro/simstack).
 %  Copyright (c) 2020 Henning Hollermann.
 %  
 %  This program is free software: you can redistribute it and/or modify  
 %  it under the terms of the GNU General Public License as published by  
 %  the Free Software Foundation, version 3.
 % 
 %  This program is distributed in the hope that it will be useful, but 
 %  WITHOUT ANY WARRANTY; without even the implied warranty of 
 %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 %  General Public License for more details.
 % 
 %  You should have received a copy of the GNU General Public License 
 %  along with this program. If not, see <http://www.gnu.org/licenses/>.
 %
 %  Example for the usage of buildStack and SimStack simulating diffraction
 %  at superlattice of 10x[6 nm Sb2Te3 / 1 nm GeTe].
 %}


clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load classes
addpath(genpath('D:\Science Data\Projects\Software Development\Matlab'));

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation range, x axis in A^-1
Qz = 1:0.0002:6;

% lattice constants of the constituents
c_Sb2Te3 = 30.5287;
c_GeTe = 10.6917;
% thicknesses of the building blocks
t_QL = c_Sb2Te3/3;
t_BL = c_GeTe/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stack definition

% number of superposed random structures
n_struct = 20;
% CSL repetitions
repetitions = 10;
% sample definition
n_QL = 6;
n_BL = 3;
% sigma about which the number of building blocks are allowed to vary
s_QL = 1.5;
s_BL = 1.5;
% calculate SL BiLayer Thickness
lambda_sim = n_QL*t_QL + n_BL*t_BL;
eta_sim = n_QL*t_QL/lambda_sim;
t_sim = lambda_sim*repetitions;

fprintf('Running %s.m: SimStack Minimal Example\n',mfilename)
fprintf(['Simulated repetition length: %g nm\nSb2Te3-content eta: %g\n'...
    'Total thickness: %g nm\n'],lambda_sim/10,eta_sim,t_sim/10)

%% Sb2Te3 QT disordered structure simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % calculate SL BiLayer Thickness
lambda_sim = n_QL*t_QL + n_BL*t_BL;
eta_sim = n_QL*3/n_BL;
t_sim = lambda_sim*repetitions;
c_GeTe=t_BL*3;
c_Sb2Te3=t_QL*3;

fprintf(['Simulated repetition length Lambda: %g nm\n', ...
'Sample: %ix[%0.2f nm GeTe / %.2f nm Sb2Te3]', ...
'eta n_QL*3/n_BL: %g\n',...
'Total thickness: %g nm\n'],...
lambda_sim/10,repetitions,n_BL*t_BL/10,n_QL*t_QL/10,eta_sim,t_sim/10)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate randomized structures with the "buildStack" class
tic
fprintf('Generating %i structures with some randomness...\n', n_struct)
for k=1:n_struct
    n = repetitions;

    stack = buildStack;
    block_BL = stack.BL(c_GeTe);
    block_QL = stack.QT_Kooi(c_Sb2Te3);
    
    % build CSL structure
    for i=1:n
        % Sb2Te3 quintuples
        n_QL_actual = -1;
        while n_QL_actual < 0; n_QL_actual = n_QL+round((randn()*sqrt(s_QL))); end
        stack.multiaddBlock(block_QL,n_QL_actual)
        % GeTe bilayers
        n_BL_actual = -1;
        while n_BL_actual < 0; n_BL_actual = n_BL+round((randn()*sqrt(s_BL))); end
        stack.multiaddBlock(block_BL,n_BL_actual)
    end
    
    % add top layer Sb2Te3
    n_QL_actual = -1;
    while n_QL_actual < 0; n_QL_actual = n_QL+round((randn()*sqrt(s_QL))); end
    stack.multiaddBlock(block_QL,n_QL_actual)
    stack.addAtom('Te',0.000);
    stack_list{k} = stack;
end
toc
% create SimStack object
SimStackQTDisordered = SimStack('qz',Qz,'corr','LPTF');
% add all the structures
for k=1:n_struct
    SimStackQTDisordered.addStructure(stack_list{k}.elements,stack_list{k}.zpos,stack_list{k}.z);
end
% calculate
SimStackQTDisordered.update;
% plot
SimStackQTDisordered.plot;