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
 %  at a perfect superlattice of 10x[8 GaAs / 4 AlAs].
 %}

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load classes SimStack and buildStack
addpath(genpath('/path/to/SimStack'))

c_GaAs = 5.6537; % in anstroms
c_AlAs = 5.6620; % in anstroms

% build perfect SL stack with 10x[8 GaAs / 4 AlAs]
SLStack = buildStack;
% used building blocks are defined as functions within the class
% they define the 'zpos' of the 'elements' in absolute coordinates
% scaled by the given c lattice parameter
% 'space' is the absolute distance to the subsequent block
block_GaAs = SLStack.GaAs(c_GaAs);
block_AlAs = SLStack.AlAs(c_AlAs);
% add GaAs 8 GaAs blocks with in total 16 atoms
SLStack.multiaddBlock(block_GaAs,8)
% add AlAs layer
SLStack.multiaddBlock(block_AlAs,4);
% simulation space
Qz=3:0.0001:6;
% create SimStack object with 10 unit cells
SLSim = SimStack('qz',Qz,'corr','none','n_UC',10);
% add SL stack
SLSim.addStructure(SLStack);
% calculate all quantities including the added structure
SLSim.update; 
% plot result
SLSim.plot;

