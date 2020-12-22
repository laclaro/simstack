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
 %  SimStack uses kinematical diffraction theory to simulate a diffraction 
 %  experiment on the Qz axis. For this purpose, the scattered intensity 
 %  from one or more chains of atoms are calculated and averaged.
 %  The chains of atoms can be setup using the buildStack class and may contain
 %  crystal defects.
 %
 %  For a background see the dissertations of Felix Lange and Henning 
 %  Hollermann at RWTH Aachen University, available at
 %  https://publications.rwth-aachen.de
 %}

classdef SimStack < handle
    properties
        Qz = []; % simulation range
        theta = [];
        lambda = 1.54056;
        zpos = {};     % cell of atom position arrays
        elements = {}; % cell of cells of element-strings
        species = {};  % elemental species
        stoichiometry = {}; % 
        F = {}; % simulated structure factors
        G = {}; % simulated crystal shape function
        I = {}; % simulated intensities
        I_coh = []; % calculated coherent sum
        I_incoh = []; % calculated incoherent sum
        corr = 'LPTF';
    end
    
    properties (Access = private)
        % periodic = false; 
        n_struct = 0;
        n_UC = 1; % number of unit cells in each stack
        size_stack = []; % stack sizes provided by the given stack objects
        struct_list = {}; % legacy
    end
    
    methods (Access = public)
   
        % setup the SimStack object, calculating all values
        function obj = SimStack(varargin)


            %  should we use struct_list?
            legacy = false;
            ind_legacy = find(strcmp(varargin,'legacy'),1);
            if ~isempty(ind_legacy)
                legacy = true;
            end
            ind_modern = find(strcmp(varargin,'modern'),1);
            if ~isempty(ind_modern)
                legacy = false;
            end
            
            % by default directly use a cell of double arrays for the
            % z-positions and a second 
            if legacy == false
                % warning('Modern SimStack is used.')
                obj.SimStackInit(varargin{:});
            % use legacy option to give a structure list as a cell of
            % buildStructure objects
            else
                rest_varargin = {};
                warning('Legacy is true, expecting cell array of buildStructure objects.')
                if length(varargin) > 1
                    rest_varargin = {varargin{2:end}};
                end
                obj.SimStackInitLegacy(varargin{1},rest_varargin{:});
            end

        end
        

        % calculate and set F, G and I for each structure in struct_list
        function simulate(obj,varargin)
            
            % default
            parpool = true;
            n_struct = length(obj.zpos);
            
            % default: simulate all
            k_simulate = 1:n_struct;
            
            % test options parpool and which structures to simulate
            i=1;
            while i <= length(varargin)
                if(ischar(varargin{i}) == 1)
                switch varargin{i}
                    case 'parpool'
                        parpool = varargin{i+1};
                        i=i+1;
                    case { 'k', 'only'}
                        k_simulate = varargin{i+1};
                        i = i+1;
                end
                end
                i = i+1;
            end
            
            % evaluate corrections to apply
            str_corr_factor = '1';
            if ~isempty(strfind(obj.corr ,'L'))
                str_corr_factor = [str_corr_factor '.*obj.calcL'];
            end
            if ~isempty(strfind(obj.corr ,'P'))
                str_corr_factor = [str_corr_factor '.*obj.calcP'];
            end
            if ~isempty(strfind(obj.corr ,'T')) || ~isempty(strfind(obj.corr ,'F'))
                str_corr_factor = [str_corr_factor '.*obj.calcTF'];
            end
            if ~isempty(strfind(obj.corr,'LPTF'))
                str_corr_factor = 'obj.calcL.*obj.calcP.*obj.calcTF';
            end
            if strcmpi(obj.corr,'none')
                corr_factor = 1;
            else
                str_corr_factor = regexprep(str_corr_factor,'^1\.\*','');
                corr_factor = eval(str_corr_factor);
            end
            
            % eventually use parpool
            if (parpool==true && length(k_simulate) > 1)
                % initialize multicore processing pool
                if isempty(gcp('nocreate')) parpool; end                
            end
            % calculate F, G and I for each structure
            fprintf('Calculate F_cell, G_cryst and I with %s corrections',obj.corr);
            local_G = repmat({''},n_struct,1);
            local_F = local_G; local_I = local_G;
            % with parfor of without, speeds up the process significantly
            if parpool == true
                fprintf(' using parpool:');
                fprintf(['\n' repmat('.',1,n_struct) '\n\n']);
%                 L = obj.calcL; P = obj.calcP; TF = obj.calcTF;
                parfor k = k_simulate
                    local_F{k} =  obj.calcF(k);
                    local_G{k} =  obj.calcG(k);
                    local_I{k} =  abs(local_F{k}.*local_G{k}.*corr_factor).^2;
                    fprintf('\b|\n');
                end
            else
                if n_struct > 1; fprintf(':\nNote: parpool may speed up the process'); end
                fprintf(['\n' repmat('.',1,n_struct) '\n\n']);
                for k = k_simulate
                    local_F{k} =  obj.calcF(k);
                    local_G{k} =  obj.calcG(k);
                    local_I{k} =  abs(local_F{k}.*local_G{k}.*corr_factor).^2;
                    fprintf('\b|\n');
                end
            end
            
            % set properties
            for k = k_simulate
                obj.F{k} = local_F{k};
                obj.G{k} = local_G{k};
                obj.I{k} = local_I{k};
            end
            obj.I_coh = obj.coh_superpos(corr_factor);
            obj.I_incoh = obj.incoh_superpos;
        end

        % get the value of I_incoh for a given Qz
        function I_incoh_val = getI_incoh(obj,Qz_arg)
            I_incoh_val = zeros(1,length(Qz_arg));
            for i=1:length(Qz_arg)
                [~, idx] = min(abs(obj.Qz-Qz_arg(i)));
                I_incoh_val(i) = obj.I_incoh(idx);
            end
        end

        % get the value of I_coh for a given Qz
        function I_coh_val = getI_coh(obj,Qz_arg)
            [~, idx] = min(abs(obj.Qz-Qz_arg));
            I_coh_val = obj.I_coh(idx);
        end
        
        % gets unique elements within the given structure of index k
        function species = getSpecies(obj,k)
            species = unique(obj.elements{k});
        end
        
        % select/filter atoms in a structure by element and return
        % zpositions and (all the same) elements
        function [filtered_elements,filtered_zpos] = filter_by_element(~,elements,zpos,filter)
            filtered_ind = ismember(elements,filter);
            % should be all the same, just to check
            filtered_elements = elements(filtered_ind);
            % zpositions of filtered elements
            filtered_zpos = zpos(filtered_ind);
        end
        
        % calculates and returns the stoichiometry within the given structure
        % k: struct_list index
        function [element, content] = getStoichiometry(obj,k)
            % number of atoms in structure
            n_atoms = length(obj.zpos{k});
            % unique elements in structure
            unique_el = obj.getSpecies(k);
            % stoichiometry
            content = zeros(1,length(unique_el));
            % number of atoms per specie in structure
            for i=1:length(unique_el)
                [atoms_el, ~] = obj.filter_by_element(obj.elements{k},obj.zpos{k},unique_el{i});
                n_atoms_el = length(atoms_el);
                content(i) = n_atoms_el/n_atoms*100;
            end
            element = unique_el;
        end
        
        % find zpositions of vdW gaps
        function zpos_vdW = find_vdW_gaps(obj,k)
            vdW_gap = 2.8;
            atom_diff = diff(obj.zpos{k});
            idx = find(atom_diff>vdW_gap);
            zpos_vdW = obj.zpos{k}(idx)+(atom_diff(idx)./2); 
        end
        
        % plot simulation data and structure
        function [fig_struct, fig_sim] = plot(obj,varargin)
            fig_struct = obj.plot_struct(varargin{:});
            if isempty(obj.I) == 0
                fig_sim = obj.plot_sim(varargin{:});
            else
                fig_sim = 0;
            end
        end
        
        % plot simulation data
        % title, 'titlestring' is the optional title
        function fig_sim = plot_sim(obj,varargin)
            
%             figure;
            cmap_coh = summer;
            cmap_incoh = winter;
%             figs = findobj('type','figure');
%             close(figs(1));
            
            hold_on = false;
            update = false;
            % by default, don't print a title
            titlestr= '';
            i=1;
            while i <= length(varargin)
                if(ischar(varargin{i}) == 1)
                switch varargin{i}
                    case 'title'
                        titlestr = varargin{i+1};
                        i=i+1;
                    case 'notitle'
                        titlestr = '';
                        i=i+1;
                    case 'update'
                        update = true;
                    case 'hold'
                        update = true;
                        hold_on = true;
                end
                end
                i = i+1;
            end
            
            % reuse the figure if update is true 
            if (update == true)
                fig_sim = findobj('type','figure','name','SimStack Simulation');
            end
            % open a new figure if "fig" is empty
            if ~exist('fig_sim','var') || isempty(fig_sim);
                fig_sim = figure('name','SimStack Simulation');
            end
            % select figure
            figure(fig_sim);

            % set the colors
            num_lines = numel(findobj(fig_sim,'Type','line'));
            color_coh = cmap_coh(mod(ceil(num_lines/2),length(cmap_coh))+1,:);
            color_incoh = cmap_incoh(mod(ceil(num_lines/2),length(cmap_incoh))+1,:);

            % add plots to figure
            if hold_on == true; hold on; else hold off; end
            
            % plot first curve
            hcoh = plot(obj.Qz,obj.I_coh,'color',color_coh);
            xlabel(['Q_z / ' char(197) '^{-1}'])
            ylabel('Intensity')
            set(gca,'YScale','log')
            if ~isempty(titlestr)
                title(titlestr);
            end
            grid on
            box on
            
            if (length(obj.zpos)) > 1
                % add second curve
                hold on
                hincoh = plot(obj.Qz,obj.I_incoh,'color',color_incoh);
                l = legend([hcoh hincoh], {...
                    sprintf('coherent superposition'), ...
                    sprintf('incoherent superposition')
                    }, 'Location', 'Northwest');
            else
                l = legend(hcoh, {...
                    sprintf('simulated intensity')
                    }, 'Location', 'Northwest');
            end
            
             hold off
        end
        
        % plot all atom columns eventually with vdW gaps
        % triggered by options ('plot', true, 'plot_vdW', true)
        function fig_struct = plot_struct(obj, varargin)

            update = false;
            % default colormap
            Colormap = @ElementColors;
            % plot all structs by default
            n_struct = length(obj.zpos);
            k_list = 1:n_struct;
            % vdW plot style
            vdW_plotstyle = {'square','color',Colormap('vdW',7)};
            % atoms plot style
            atoms_plotstyle = {'MarkerSize',3};
            % plot vdW gaps?
            plot_vdW = false;
            i=1;
            while i <= length(varargin)
                if(ischar(varargin{i}) == 1)
                switch varargin{i}
                    case {'element_colors'}
                        Colormap = varargin{i+1};
                        i=i+1;
                    case 'k'
                        k_list = varargin{i+1};
                        i=i+1;
                    case {'vdW', 'plot_vdW'}
                        plot_vdW = varargin{i+1};
                        i=i+1;
                    case {'update','hold'}
                        update = true;
                end
                end
                i = i+1;
            end
            
            % reuse the figure if update is true
            if (update == true)
                fig_struct = findobj('type','figure','name','SimStack Structures');
            end
            % open a new figure if "fig" is empty
            if ~exist('fig_struct','var') || isempty(fig_struct);
                fig_struct = figure('name','SimStack Structures');
            end
            % select figure
            figure(fig_struct);

            xlabel('atom column number n');
            ylabel(['z / ' char(197)]);
            % find atoms
            hold on;
            for i=numel(obj.species):-1:1
                zpos = []; struct_pos = [];
                ele = {};
                for k = k_list
                    [~, tmp_zpos] = obj.filter_by_element(obj.elements{k},obj.zpos{k},obj.species{i});
                    struct_pos = [struct_pos; k*ones(length(tmp_zpos),1)];
                    zpos = [zpos; tmp_zpos];
                end
                % in case, plot vdW gaps below
                if plot_vdW == true
                    zpos_gaps = [];
                    struct_pos_gaps = [];
                    for k = 1:n_struct
                        tmp_zpos_gaps = obj.find_vdW_gaps(k);
                        zpos_gaps = [zpos_gaps; tmp_zpos_gaps];
                        struct_pos_gaps = [struct_pos_gaps; k*ones(length(tmp_zpos_gaps),1)];
                    end
                    plot(struct_pos_gaps,zpos_gaps,vdW_plotstyle{:})
                end
                % plot the atoms
                plot(struct_pos,zpos,'o',...
                    'MarkerEdgeColor',Colormap(obj.species{i},3),...
                    'MarkerFaceColor',Colormap(obj.species{i},1),...
                    atoms_plotstyle{:});
            end
            % setup the legend
            legend_handles = [];
            legend_texts = {''};
            for i=1:numel(obj.species)
                legend_handles(i) = plot(0,0,'o',...
                    'MarkerEdgeColor',Colormap(obj.species{i},3),...
                    'MarkerFaceColor',Colormap(obj.species{i},1),...
                    atoms_plotstyle{:},'visible','off','DisplayName',obj.species{i});
                legend_texts{i} = obj.species{i};
            end
            % add vdW Gaps to legend
            if plot_vdW == true
                legend_handles(end+1) = plot(0,0,vdW_plotstyle{:},'visible','off');
                legend_texts{end+1} = 'vdW';
            end
            legend(legend_handles,legend_texts);
            
            hold off
        end
        
        % write simulated data to file
        function export_data(obj,filename)
            
            A = [obj.Qz; obj.I_coh; obj.I_incoh]';
            
            [dir,mfile,~] = fileparts(mfilename('fullpath'));
            
            fid = fopen(filename,'wt');
            % write header
            header = sprintf('# File created by SimStack (%s) on %s\n# %s;%s;%s\n',...
                datestr(datetime), mfile,'Qz','I_coh','I_incoh');
            fprintf(fid, header);
            % have a matrix A ready to go, append data
            dlmwrite(filename, A,'delimiter', ';', '-append')
            fclose(fid);
        end

        % add structure to the list
        function addStructure(obj,varargin)
            if isa(varargin{1},'buildStack')
                % add another atom coloumn
                obj.elements{end+1} = varargin{1}.elements';
                obj.zpos{end+1} = varargin{1}.zpos';
                obj.size_stack(end+1) = varargin{1}.z;
            else
                % add another atom coloumn
                obj.elements{end+1} = varargin{1}'; % elements
                obj.zpos{end+1} = varargin{2}'; % zpos
                obj.size_stack(end+1) = varargin{3}; % stack_size / z
            end
            
            if ~isempty(find(strcmp(varargin,'update'),1))
                % update because structure was added
                obj.update();
            end
        end
        
        % add structure to the list
        function delStructure(obj,k,varargin)

            % delete another atom coloumn
            obj.elements(k) = [];
            obj.zpos(k) = [];
            
            if ~isempty(find(strcmp(varargin,'update'),1))
                % update because structure was added
                obj.update();
            end
        end
        
        % update variables whenever a new structure is added
        function update(obj)
            
            n_struct = length(obj.zpos);

            for k = 1:n_struct
                % get unique elements
                obj.species = sort(union(obj.species,obj.getSpecies(k)));
                % calculate per-struct-stoichiometries
                [obj.stoichiometry{k}{1}, obj.stoichiometry{k}{2}] = obj.getStoichiometry(k);
            end
            
            % calculate values
            if length(obj.F) == (n_struct-1)
                obj.simulate('parpool', false, 'k', n_struct);
            else
                obj.simulate('parpool', true);
            end
            
        end
        
        % smooth y data by using a movin average
        function smooth(obj, amount)
%             n = length(obj.y);
%             % FFT
%             f = fft(obj.y);
%             % cut some frequencies
%             f(n/2+1-amount:n/2+amount) = zeros(2*amount,1);
%             % inverse FFT
%             obj.y = real(ifft(f));
            obj.I_coh = smooth(obj.I_coh, amount, 'moving');
            obj.I_incoh = smooth(obj.I_incoh, amount, 'moving');
            
        end
        
        % calculate the structure factor ensemble average
        % atompos is generated by the buildStructure class
        %{
        function F_cell_av = calcSteinerI(obj,lambda)
            
            for k=1:obj.n_struct
                % calculate atomic form factors for unique elements in structure
                unique_el = obj.getSpecies(k);
                unique_el_fs = ones(length(unique_el),length(obj.Qz));
                for i=1:length(unique_el)
                    unique_el_fs(i,:) = obj.atomic_form_factor(unique_el(i),obj.Qz);
                end

                % initialize and calculate the structure factor F_cell
                % parfor does not speed up the process here
                F_cell = zeros(1,length(obj.Qz));
                for s=1:length(obj.zpos{k})
                    % match element name without trailing number
                    match = regexp(obj.elements{k}{s},'[a-zA-Z]+','match');
                    % make all atoms Tellurium
                    %match = regexp('Te','[a-zA-Z]+','match');
                    % index of the element to select form factor (lower for "non-case-sensitivity")
                    ind = ismember(lower(unique_el),lower(match{1}));
                    fs = unique_el_fs(ind,:);

                    % only in qz direction, alpha
                    F_cell = F_cell + fs.*exp(-1i*(obj.Qz.*obj.zpos{k}(s)));
                    alpha_norm = alpha_norm + exp(-1i*(obj.Qz.*obj.zpos{k}(s)));
                    
                    gamma = gamma + abs(fs.*exp(-1i*(obj.Qz.*obj.zpos{k}(s)))).^2;
                    gamma_norm = gamma_norm + abs(exp(-1i*(obj.Qz.*obj.zpos{k}(s)))).^2;
                    
                    xi = xi + exp(i.*Qz.*lambda);
                    xi_norm = xi_norm + exp(i.*Qz.*lambda);
                    
    %                 F_cell = F_cell + fs.*exp(-1i*(Qz.*mod(zpositions(s)+UC_shift*size_UC,size_UC)));
                end
                
                % beta, conj F_cell
                beta = beta + conj(F_cell) .* exp(-i*obj.Qz*lambda);
                
                alpha_sum = alpha_sum + F_cell;
                alpha_norm_sum = alpha_norm_sum + alpha_norm;
                
                beta_sum = 
                beta_norm = obj.n_struct * exp(-i*obj.Qz*lambda);
                gamma = 
            end
            
            % structure factor ensemble average
            alpha = alpha_sum./alpha_norm_sum;
        end
        %}
        
        
    end

    
    
    % methods that are not accessible directly (helper-functions)
    methods (Access = private)
    
        function SimStackInit(obj,varargin)
            
            % DEFAULTS
            simulate = true;
            plot = false;
            plot_vdW = false;
            parpool = true;
            export_data = false;
            export_filename='SimStack_data.csv';
            % default simulation range
            obj.Qz = 1:0.001:5;

            
            i=1;
            while i <= length(varargin)
                if(ischar(varargin{i}) == 1)
                switch varargin{i}
                    case 'parpool'
                        parpool = varargin{i+1};
                        i=i+1;
                    case {'n', 'n_UC'}
                        obj.n_UC = varargin{i+1};
                        i=i+1;
                    case 'plot'
                        plot = varargin{i+1};
                        i = i+1;
                    case 'simulate'
                        simulate = varargin{i+1};
                        i = i+1;
                    case 'export_data'
                        export_data = varargin{i+1};
                        i = i+1;
                    case 'filename'
                        export_filename = varargin{i+1};
                        i = i+1;
                    case 'qz'
                        obj.Qz = varargin{i+1};
                        i = i+1;
                    case {'zpos', 'pos'}
                        obj.elements = varargin{i+1};
                        i = i+1;
                    case {'atoms', 'elements'}
                        obj.zpos = varargin{i+1};
                        i = i+1;
                    case 'corr'
                        obj.corr = varargin{i+1};
                        i = i+1;
                end
                end
                i = i+1;
            end
            
            obj.theta = asin(obj.Qz.*obj.lambda./(4*pi));
            
%              % number of structures
%              obj.n_struct = length(obj.elements);
                
%             % assign n_UC
%             if isempty(obj.n_UC)
%               obj.n_UC = 1;
%             end
            
        end

       
        % expect buildStructure as struc_list
        function SimStackInitLegacy(obj,struct_list,varargin)
  
            % DEFAULTS
            simulate = true;
            plot = false;
            plot_vdW = false;
            parpool = true;
            export_data = false;
            export_filename='SimStack_data.csv';
%             % default simulation range
            obj.Qz = 1:0.001:5;
            % number of structures
            obj.n_struct = length(struct_list);
            
            % test options
            i=1;
            while i <= length(varargin)
                if(ischar(varargin{i}) == 1)
                switch varargin{i}
                    case 'parpool'
                        parpool = varargin{i+1};
                        i=i+1;
%                     case 'periodic'
%                         obj.periodic = varargin{i+1};
%                         i=i+1;
                    case 'n'
                        obj.n_UC = varargin{i+1};
                        i=i+1;
                    case 'plot'
                        plot = varargin{i+1};
                        i = i+1;
                    case 'simulate'
                        simulate = varargin{i+1};
                        i = i+1;
                    case 'export_data'
                        export_data = varargin{i+1};
                        i = i+1;
                    case 'filename'
                        export_filename = varargin{i+1};
                        i = i+1;
                    case {'XLim', 'qz', 'range'}
                        obj.Qz = varargin{i+1};
                        i = i+1;
                end
                end
                i = i+1;
            end
            
            % check if there was given a single structure or a list of
            % structures
            if isa(struct_list,'buildStructure')
                disp('struct_list is a single buildStructure object')
                obj.struct_list = {struct_list};
                obj.convert_struct_list;
            elseif isa(struct_list,'buildStack')
                disp('struct_list is a single buildStack object')
                obj.struct_list = {struct_list};
                obj.convert_stack_list;
            elseif isa(struct_list,'cell')
                if isa(struct_list{1},'buildStructure')
                    disp('struct_list is a cell of buildStructure objects')
                    obj.struct_list = struct_list;
                    obj.convert_struct_list;
                elseif isa(struct_list{1},'buildStack')
                    disp('struct_list is a cell of buildStack objects')
                    obj.struct_list = struct_list;
                    obj.convert_stack_list
                else
                    error('ERROR: given structure list ist incompatible with SimStack')
                end
            else
                error('ERROR: First argument must be a buildStructure, buildStack or a cell of such!')
            end

%             % sample thickness
%             if obj.periodic == false
%                 obj.n_UC = 1;
%             end
            
            for k = 1:obj.n_struct
                % get unique elements
                obj.species = sort(union(obj.species,obj.getSpecies(k)));
                % calculate per-struct-stoichiometries
                [obj.stoichiometry{k}{1}, obj.stoichiometry{k}{2}] = obj.getStoichiometry(k);
            end
            
            % calculate values
            if simulate == true
                obj.simulate('parpool', parpool);
            else
                fprintf(2,'WARNING: Since option "simulate" is set to false, diffraction intensities\n are not calculated.\n');
            end
            
            % plot structure
            if plot == true
                [fig_struct, fig_sim] = obj.plot(varargin{:});
            end
            
            % export data
            if simulate == true && export_data == true
                obj.export_data(export_filename);
            end
            
        end

        
        
        % calculate the coherent superposition of the structures
        function I_coh_sum = coh_superpos(obj,corr_factor)
            
            n_struct = length(obj.zpos);
            coh_sum = zeros(1,length(obj.Qz));
            for k = 1:n_struct
                coh_sum = coh_sum + obj.G{k}.*obj.F{k};
            end
            I_coh_sum = abs(coh_sum.*corr_factor).^2;
            obj.I_coh = I_coh_sum;
        end
        
        % calculate the incoherent superpositions of all atom columns  
        function I_incoh_sum = incoh_superpos(obj)
            I_incoh_sum = sum(cell2mat(obj.I'),1);
            obj.I_incoh = I_incoh_sum;
        end
        
        % calculate the structure factor of the unit cell of kth structure
        % atompos is generated by the buildStructure class
        function F_cell = calcF(obj,k)

            % calculate atomic form factors for unique elements in structure
            unique_el = obj.getSpecies(k);
            unique_el_fs = ones(length(unique_el),length(obj.Qz));
            for i=1:length(unique_el)
                unique_el_fs(i,:) = obj.atomic_form_factor(unique_el(i),obj.Qz);
            end

            % initialize and calculate the structure factor F_cell
            % parfor does not speed up the process here
            F_cell = zeros(1,length(obj.Qz));
            for s=1:length(obj.zpos{k})
                % match element name without trailing number
                match = regexp(obj.elements{k}{s},'[a-zA-Z]+','match');
                % make all atoms Tellurium
                %match = regexp('Te','[a-zA-Z]+','match');
                % index of the element to select form factor (lower for "non-case-sensitivity")
                ind = ismember(lower(unique_el),lower(match{1}));
                fs = unique_el_fs(ind,:);

                % only in qz direction
                F_cell = F_cell + fs.*exp(-1i*(obj.Qz.*obj.zpos{k}(s)));
%                 F_cell = F_cell + fs.*exp(-1i*(Qz.*mod(zpositions(s)+UC_shift*size_UC,size_UC)));
            end
        end
        

        % calculate laue function for kth structure
        function G_cryst = calcG(obj,k)
%             size_UC = obj.size_stack(k)
%             G_cryst = exp(1i*obj.Qz*size_UC/2).*sin(size_UC*obj.n_UC*obj.Qz/2)./sin(size_UC*obj.Qz/2);
            G_cryst = (exp(-1i.*obj.Qz.*obj.n_UC.*obj.size_stack(k))-1)./(exp(-1i.*obj.Qz.*obj.size_stack(k))-1);
%             exp(-1i.*obj.Qz.*obj.size_stack(k)).*
        end
        
        % calculate Lorentz Polarization factor as described in Pecharsky's
        % "FUNDAMENTALS OF POWDER DIFFRACTION AND STRUCTURAL CHARACTERIZATION 
        % OF MATERIALS" in chapter 2.10.4
        function L = calcL(obj)
            L = 1./(4*cos(obj.theta).*sin(obj.theta).^2);
        end
        
        function P = calcP(obj)
            P = (1+cos(2.*obj.theta).^2)/2;
        end

        
        % calculate Temperature factor as described in Pecharsky's
        % "FUNDAMENTALS OF POWDER DIFFRACTION AND STRUCTURAL CHARACTERIZATION 
        % OF MATERIALS" in chapter 2.11.3 
        function TF = calcTF(obj)
            B=1;
            TF = exp(-B*sin(obj.theta).^2/obj.lambda);
        end
        
        % getting approximated curve for atomic form factors
        % atomic_form_factor.m should be only run once per unique element in the
        % structure and applied for every atom of that type in the structure.
        function f = atomic_form_factor(~,element,q)

            fstr = ['c(1).*exp(-c(2).*(qsquare./(4*pi).^2)) +' ...
                'c(3).*exp(-c(4).*(qsquare./(4*pi).^2))  +' ...
                'c(5).*exp(-c(6).*(qsquare./(4*pi).^2))  +' ...
                'c(7).*exp(-c(8).*(qsquare./(4*pi).^2))  +' ...
                'c(9)'];

            % do a proper dot product in case q is a vector
            qsquare=q.*q;

            % element names must coincide with c_values from
            % http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

            % no significant speed boost if reduced list is used
            % c_elements = {'Ge', 'Sb', 'Te','Ga','As'};
            % 
            % c_values = {[16.0816 2.8509 6.3747 0.2516 3.7068 11.4468 3.683 54.7625 2.1313], ...
            %      [19.6418 5.3034 19.0455 0.4607 5.0371 27.9074 2.6827 75.2825 4.5909], ...
            %      [19.9644 4.81742 19.0138 0.420885 6.14487 28.5284 2.5239 70.8403 4.352], ...
            %      [15.2354 3.0669 6.7006 0.2412 4.3591 10.7805 2.9623 61.4135 1.7189], ...
            %      [16.6723 2.6345 6.0701 0.2647 3.4313 12.9479 4.2779 47.7972 2.531]};

            c_elements = {'H','H1-','He','Li','Li1+','Be','Be2+','B','C','Cval','N',...
                'O','O1-','F','F1-','Ne','Na','Na1+','Mg','Mg2+','Al','Al3+',...
                'Siv','Sival','Si4+','P','S','Cl','Cl1-','Ar','K','K1+','Ca','Ca2+',...
                'Sc','Sc3+','Ti','Ti2+','Ti3+','Ti4+','V','V2+','V3+','V5+','Cr',...
                'Cr2+','Cr3+','Mn','Mn2+','Mn3+','Mn4+','Fe','Fe2+','Fe3+','Co',...
                'Co2+','Co3+','Ni','Ni2+','Ni3+','Cu','Cu1+','Cu2+','Zn','Zn2+',...
                'Ga','Ga3+','Ge','Ge4+','As','Se','Br','Br1-','Kr','Rb','Rb1+','Sr',...
                'Sr2+','Y','Y3+','Zr','Zr4+','Nb','Nb3+','Nb5+','Mo','Mo3+','Mo5+',...
                'Mo6+','Tc','Ru','Ru3+','Ru4+','Rh','Rh3+','Rh4+','Pd','Pd2+','Pd4+',...
                'Ag','Ag1+','Ag2+','Cd','Cd2+','In','In3+','Sn','Sn2+','Sn4+','Sb',...
                'Sb3+','Sb5+','Te','I','I1-','Xe','Cs','Cs1+','Ba','Ba2+','La',...
                'La3+','Ce','Ce3+','Ce4+','Pr','Pr3+','Pr4+','Nd','Nd3+','Pm',...
                'Pm3+','Sm','Sm3+','Eu','Eu2+','Eu3+','Gd','Gd3+','Tb','Tb3+','Dy',...
                'Dy3+','Ho','Ho3+','Er','Er3+','Tm','Tm3+','Yb','Yb2+','Yb3+','Lu',...
                'Lu3+','Hf','Hf4+','Ta','Ta5+','W','W6+','Re','Os','Os4+','Ir',...
                'Ir3+','Ir4+','Pt','Pt2+','Pt4+','Au','Au1+','Au3+','Hg','Hg1+',...
                'Hg2+','Tl','Tl1+','Tl3+','Pb','Pb2+','Pb4+','Bi','Bi3+','Bi5+',...
                'Po','At','Rn','Fr','Ra','Ra2+','Ac','Ac3+','Th','Th4+','Pa','U',...
                'U3+','U4+','U6+','Np','Np3+','Np4+','Np6+','Pu','Pu3+','Pu4+',...
                'Pu6+','Am','Cm','Bk','Cf','vacuum'};
            c_values = {[0.489918	20.6593	0.262003	7.74039	0.196767	49.5519	0.049879	2.20159	0.001305], ...
            [0.897661	53.1368	0.565616	15.187	0.415815	186.576	0.116973	3.56709	0.002389], ...
            [0.8734	9.1037	0.6309	3.3568	0.3112	22.9276	0.178	0.9821	0.0064], ...
            [1.1282	3.9546	0.7508	1.0524	0.6175	85.3905	0.4653	168.261	0.0377], ...
            [0.6968	4.6237	0.7888	1.9557	0.3414	0.6316	0.1563	10.0953	0.0167], ...
            [1.5919	43.6427	1.1278	1.8623	0.5391	103.483	0.7029	0.542	0.0385], ...
            [6.2603	0.0027	0.8849	0.8313	0.7993	2.2758	0.1647	5.1146	-6.1092], ...
            [2.0545	23.2185	1.3326	1.021	1.0979	60.3498	0.7068	0.1403	-0.1932], ...
            [2.31	20.8439	1.02	10.2075	1.5886	0.5687	0.865	51.6512	0.2156], ...
            [2.26069	22.6907	1.56165	0.656665	1.05075	9.75618	0.839259	55.5949	0.286977], ...
            [12.2126	0.0057	3.1322	9.8933	2.0125	28.9975	1.1663	0.5826	-11.529], ...
            [3.0485	13.2771	2.2868	5.7011	1.5463	0.3239	0.867	32.9089	0.2508], ...
            [4.1916	12.8573	1.63969	4.17236	1.52673	47.0179	-20.307	-0.01404	21.9412], ...
            [3.5392	10.2825	2.6412	4.2944	1.517	0.2615	1.0243	26.1476	0.2776], ...
            [3.6322	5.27756	3.51057	14.7353	1.26064	0.442258	0.940706	47.3437	0.653396], ...
            [3.9553	8.4042	3.1125	3.4262	1.4546	0.2306	1.1251	21.7184	0.3515], ...
            [4.7626	3.285	3.1736	8.8422	1.2674	0.3136	1.1128	129.424	0.676], ...
            [3.2565	2.6671	3.9362	6.1153	1.3998	0.2001	1.0032	14.039	0.404], ...
            [5.4204	2.8275	2.1735	79.2611	1.2269	0.3808	2.3073	7.1937	0.8584], ...
            [3.4988	2.1676	3.8378	4.7542	1.3284	0.185	0.8497	10.1411	0.4853], ...
            [6.4202	3.0387	1.9002	0.7426	1.5936	31.5472	1.9646	85.0886	1.1151], ...
            [4.17448	1.93816	3.3876	4.14553	1.20296	0.228753	0.528137	8.28524	0.706786], ...
            [6.2915	2.4386	3.0353	32.3337	1.9891	0.6785	1.541	81.6937	1.1407], ...
            [5.66269	2.6652	3.07164	38.6634	2.62446	0.916946	1.3932	93.5458	1.24707], ...
            [4.43918	1.64167	3.20345	3.43757	1.19453	0.2149	0.41653	6.65365	0.746297], ...
            [6.4345	1.9067	4.1791	27.157	1.78	0.526	1.4908	68.1645	1.1149], ...
            [6.9053	1.4679	5.2034	22.2151	1.4379	0.2536	1.5863	56.172	0.8669], ...
            [11.4604	0.0104	7.1964	1.1662	6.2556	18.5194	1.6455	47.7784	-9.5574], ...
            [18.2915	0.0066	7.2084	1.1717	6.5337	19.5424	2.3386	60.4486	-16.378], ...
            [7.4845	0.9072	6.7723	14.8407	0.6539	43.8983	1.6442	33.3929	1.4445], ...
            [8.2186	12.7949	7.4398	0.7748	1.0519	213.187	0.8659	41.6841	1.4228], ...
            [7.9578	12.6331	7.4917	0.7674	6.359	-0.002	1.1915	31.9128	-4.9978], ...
            [8.6266	10.4421	7.3873	0.6599	1.5899	85.7484	1.0211	178.437	1.3751], ...
            [15.6348	-0.0074	7.9518	0.6089	8.4372	10.3116	0.8537	25.9905	-14.875], ...
            [9.189	9.0213	7.3679	0.5729	1.6409	136.108	1.468	51.3531	1.3329], ...
            [13.4008	0.29854	8.0273	7.9629	1.65943	-0.28604	1.57936	16.0662	-6.6667], ...
            [9.7595	7.8508	7.3558	0.5	1.6991	35.6338	1.9021	116.105	1.2807], ...
            [9.11423	7.5243	7.62174	0.457585	2.2793	19.5361	0.087899	61.6558	0.897155], ...
            [17.7344	0.22061	8.73816	7.04716	5.25691	-0.15762	1.92134	15.9768	-14.652], ...
            [19.5114	0.178847	8.23473	6.67018	2.01341	-0.29263	1.5208	12.9464	-13.28], ...
            [10.2971	6.8657	7.3511	0.4385	2.0703	26.8938	2.0571	102.478	1.2199], ...
            [10.106	6.8818	7.3541	0.4409	2.2884	20.3004	0.0223	115.122	1.2298], ...
            [9.43141	6.39535	7.7419	0.383349	2.15343	15.1908	0.016865	63.969	0.656565], ...
            [15.6887	0.679003	8.14208	5.40135	2.03081	9.97278	-9.576	0.940464	1.7143], ...
            [10.6406	6.1038	7.3537	0.392	3.324	20.2626	1.4922	98.7399	1.1832], ...
            [9.54034	5.66078	7.7509	0.344261	3.58274	13.3075	0.509107	32.4224	0.616898], ...
            [9.6809	5.59463	7.81136	0.334393	2.87603	12.8288	0.113575	32.8761	0.518275], ...
            [11.2819	5.3409	7.3573	0.3432	3.0193	17.8674	2.2441	83.7543	1.0896], ...
            [10.8061	5.2796	7.362	0.3435	3.5268	14.343	0.2184	41.3235	1.0874], ...
            [9.84521	4.91797	7.87194	0.294393	3.56531	10.8171	0.323613	24.1281	0.393974], ...
            [9.96253	4.8485	7.97057	0.283303	2.76067	10.4852	0.054447	27.573	0.251877], ...
            [11.7695	4.7611	7.3573	0.3072	3.5222	15.3535	2.3045	76.8805	1.0369], ...
            [11.0424	4.6538	7.374	0.3053	4.1346	12.0546	0.4399	31.2809	1.0097], ...
            [11.1764	4.6147	7.3863	0.3005	3.3948	11.6729	0.0724	38.5566	0.9707], ...
            [12.2841	4.2791	7.3409	0.2784	4.0034	13.5359	2.3488	71.1692	1.0118], ...
            [11.2296	4.1231	7.3883	0.2726	4.7393	10.2443	0.7108	25.6466	0.9324], ...
            [10.338	3.90969	7.88173	0.238668	4.76795	8.35583	0.725591	18.3491	0.286667], ...
            [12.8376	3.8785	7.292	0.2565	4.4438	12.1763	2.38	66.3421	1.0341], ...
            [11.4166	3.6766	7.4005	0.2449	5.3442	8.873	0.9773	22.1626	0.8614], ...
            [10.7806	3.5477	7.75868	0.22314	5.22746	7.64468	0.847114	16.9673	0.386044], ...
            [13.338	3.5828	7.1676	0.247	5.6158	11.3966	1.6735	64.8126	1.191], ...
            [11.9475	3.3669	7.3573	0.2274	6.2455	8.6625	1.5578	25.8487	0.89], ...
            [11.8168	3.37484	7.11181	0.244078	5.78135	7.9876	1.14523	19.897	1.14431], ...
            [14.0743	3.2655	7.0318	0.2333	5.1652	10.3163	2.41	58.7097	1.3041], ...
            [11.9719	2.9946	7.3862	0.2031	6.4668	7.0826	1.394	18.0995	0.7807], ...
            [15.2354	3.0669	6.7006	0.2412	4.3591	10.7805	2.9623	61.4135	1.7189], ...
            [12.692	2.81262	6.69883	0.22789	6.06692	6.36441	1.0066	14.4122	1.53545], ...
            [16.0816	2.8509	6.3747	0.2516	3.7068	11.4468	3.683	54.7625	2.1313], ...
            [12.9172	2.53718	6.70003	0.205855	6.06791	5.47913	0.859041	11.603	1.45572], ...
            [16.6723	2.6345	6.0701	0.2647	3.4313	12.9479	4.2779	47.7972	2.531], ...
            [17.0006	2.4098	5.8196	0.2726	3.9731	15.2372	4.3543	43.8163	2.8409], ...
            [17.1789	2.1723	5.2358	16.5796	5.6377	0.2609	3.9851	41.4328	2.9557], ...
            [17.1718	2.2059	6.3338	19.3345	5.5754	0.2871	3.7272	58.1535	3.1776], ...
            [17.3555	1.9384	6.7286	16.5623	5.5493	0.2261	3.5375	39.3972	2.825], ...
            [17.1784	1.7888	9.6435	17.3151	5.1399	0.2748	1.5292	164.934	3.4873], ...
            [17.5816	1.7139	7.6598	14.7957	5.8981	0.1603	2.7817	31.2087	2.0782], ...
            [17.5663	1.5564	9.8184	14.0988	5.422	0.1664	2.6694	132.376	2.5064], ...
            [18.0874	1.4907	8.1373	12.6963	2.5654	24.5651	-34.193	-0.0138	41.4025], ...
            [17.776	1.4029	10.2946	12.8006	5.72629	0.125599	3.26588	104.354	1.91213], ...
            [17.9268	1.35417	9.1531	11.2145	1.76795	22.6599	-33.108	-0.01319	40.2602], ...
            [17.8765	1.27618	10.948	11.916	5.41732	0.117622	3.65721	87.6627	2.06929], ...
            [18.1668	1.2148	10.0562	10.1483	1.01118	21.6054	-2.6479	-0.10276	9.41454], ...
            [17.6142	1.18865	12.0144	11.766	4.04183	0.204785	3.53346	69.7957	3.75591], ...
            [19.8812	0.019175	18.0653	1.13305	11.0177	10.1621	1.94715	28.3389	-12.912], ...
            [17.9163	1.12446	13.3417	0.028781	10.799	9.28206	0.337905	25.7228	-6.3934], ...
            [3.7025	0.2772	17.2356	1.0958	12.8876	11.004	3.7429	61.6584	4.3875], ...
            [21.1664	0.014734	18.2017	1.03031	11.7423	9.53659	2.30951	26.6307	-14.421], ...
            [21.0149	0.014345	18.0992	1.02238	11.4632	8.78809	0.740625	23.3452	-14.316], ...
            [17.8871	1.03649	11.175	8.48061	6.57891	0.058881	0	0	0.344941], ...
            [19.1301	0.864132	11.0948	8.14487	4.64901	21.5707	2.71263	86.8472	5.40428], ...
            [19.2674	0.80852	12.9182	8.43467	4.86337	24.7997	1.56756	94.2928	5.37874], ...
            [18.5638	0.847329	13.2885	8.37164	9.32602	0.017662	3.00964	22.887	-3.1892], ...
            [18.5003	0.844582	13.1787	8.12534	4.71304	0.36495	2.18535	20.8504	1.42357], ...
            [19.2957	0.751536	14.3501	8.21758	4.73425	25.8749	1.28918	98.6062	5.328], ...
            [18.8785	0.764252	14.1259	7.84438	3.32515	21.2487	-6.1989	-0.01036	11.8678], ...
            [18.8545	0.760825	13.9806	7.62436	2.53464	19.3317	-5.6526	-0.0102	11.2835], ...
            [19.3319	0.698655	15.5017	7.98929	5.29537	25.2052	0.605844	76.8986	5.26593], ...
            [19.1701	0.696219	15.2096	7.55573	4.32234	22.5057	0	0	5.2916], ...
            [19.2493	0.683839	14.79	7.14833	2.89289	17.9144	-7.9492	0.005127	13.0174], ...
            [19.2808	0.6446	16.6885	7.4726	4.8045	24.6605	1.0463	99.8156	5.179], ...
            [19.1812	0.646179	15.9719	7.19123	5.27475	21.7326	0.357534	66.1147	5.21572], ...
            [19.1643	0.645643	16.2456	7.18544	4.3709	21.4072	0	0	5.21404], ...
            [19.2214	0.5946	17.6444	6.9089	4.461	24.7008	1.6029	87.4825	5.0694], ...
            [19.1514	0.597922	17.2535	6.80639	4.47128	20.2521	0	0	5.11937], ...
            [19.1624	0.5476	18.5596	6.3776	4.2948	25.8499	2.0396	92.8029	4.9391], ...
            [19.1045	0.551522	18.1108	6.3247	3.78897	17.3595	0	0	4.99635], ...
            [19.1889	5.8303	19.1005	0.5031	4.4585	26.8909	2.4663	83.9571	4.7821], ...
            [19.1094	0.5036	19.0548	5.8378	4.5648	23.3752	0.487	62.2061	4.7861], ...
            [18.9333	5.764	19.7131	0.4655	3.4182	14.0049	0.0193	-0.7583	3.9182], ...
            [19.6418	5.3034	19.0455	0.4607	5.0371	27.9074	2.6827	75.2825	4.5909], ...
            [18.9755	0.467196	18.933	5.22126	5.10789	19.5902	0.288753	55.5113	4.69626], ...
            [19.8685	5.44853	19.0302	0.467973	2.41253	14.1259	0	0	4.69263], ...
            [19.9644	4.81742	19.0138	0.420885	6.14487	28.5284	2.5239	70.8403	4.352], ...
            [20.1472	4.347	18.9949	0.3814	7.5138	27.766	2.2735	66.8776	4.0712], ...
            [20.2332	4.3579	18.997	0.3815	7.8069	29.5259	2.8868	84.9304	4.0714], ...
            [20.2933	3.9282	19.0298	0.344	8.9767	26.4659	1.99	64.2658	3.7118], ...
            [20.3892	3.569	19.1062	0.3107	10.662	24.3879	1.4953	213.904	3.3352], ...
            [20.3524	3.552	19.1278	0.3086	10.2821	23.7128	0.9615	59.4565	3.2791], ...
            [20.3361	3.216	19.297	0.2756	10.888	20.2073	2.6959	167.202	2.7731], ...
            [20.1807	3.21367	19.1136	0.28331	10.9054	20.0558	0.77634	51.746	3.02902], ...
            [20.578	2.94817	19.599	0.244475	11.3727	18.7726	3.28719	133.124	2.14678], ...
            [20.2489	2.9207	19.3763	0.250698	11.6323	17.8211	0.336048	54.9453	2.4086], ...
            [21.1671	2.81219	19.7695	0.226836	11.8513	17.6083	3.33049	127.113	1.86264], ...
            [20.8036	2.77691	19.559	0.23154	11.9369	16.5408	0.612376	43.1692	2.09013], ...
            [20.3235	2.65941	19.8186	0.21885	12.1233	15.7992	0.144583	62.2355	1.5918], ...
            [22.044	2.77393	19.6697	0.222087	12.3856	16.7669	2.82428	143.644	2.0583], ...
            [21.3727	2.6452	19.7491	0.214299	12.1329	15.323	0.97518	36.4065	1.77132], ...
            [20.9413	2.54467	20.0539	0.202481	12.4668	14.8137	0.296689	45.4643	1.24285], ...
            [22.6845	2.66248	19.6847	0.210628	12.774	15.885	2.85137	137.903	1.98486], ...
            [21.961	2.52722	19.9339	0.199237	12.12	14.1783	1.51031	30.8717	1.47588], ...
            [23.3405	2.5627	19.6095	0.202088	13.1235	15.1009	2.87516	132.721	2.02876], ...
            [22.5527	2.4174	20.1108	0.185769	12.0671	13.1275	2.07492	27.4491	1.19499], ...
            [24.0042	2.47274	19.4258	0.196451	13.4396	14.3996	2.89604	128.007	2.20963], ...
            [23.1504	2.31641	20.2599	0.174081	11.9202	12.1571	2.71488	24.8242	0.954586], ...
            [24.6274	2.3879	19.0886	0.1942	13.7603	13.7546	2.9227	123.174	2.5745], ...
            [24.0063	2.27783	19.9504	0.17353	11.8034	11.6096	3.87243	26.5156	1.36389], ...
            [23.7497	2.22258	20.3745	0.16394	11.8509	11.311	3.26503	22.9966	0.759344], ...
            [25.0709	2.25341	19.0798	0.181951	13.8518	12.9331	3.54545	101.398	2.4196], ...
            [24.3466	2.13553	20.4208	0.155525	11.8708	10.5782	3.7149	21.7029	0.645089], ...
            [25.8976	2.24256	18.2185	0.196143	14.3167	12.6648	2.95354	115.362	3.58324], ...
            [24.9559	2.05601	20.3271	0.149525	12.2471	10.0499	3.773	21.2773	0.691967], ...
            [26.507	2.1802	17.6383	0.202172	14.5596	12.1899	2.96577	111.874	4.29728], ...
            [25.5395	1.9804	20.2861	0.143384	11.9812	9.34972	4.50073	19.581	0.68969], ...
            [26.9049	2.07051	17.294	0.19794	14.5583	11.4407	3.63837	92.6566	4.56796], ...
            [26.1296	1.91072	20.0994	0.139358	11.9788	8.80018	4.93676	18.5908	0.852795], ...
            [27.6563	2.07356	16.4285	0.223545	14.9779	11.3604	2.98233	105.703	5.92046], ...
            [26.722	1.84659	19.7748	0.13729	12.1506	8.36225	5.17379	17.8974	1.17613], ...
            [28.1819	2.02859	15.8851	0.238849	15.1542	10.9975	2.98706	102.961	6.75621], ...
            [27.3083	1.78711	19.332	0.136974	12.3339	7.96778	5.38348	17.2922	1.63929], ...
            [28.6641	1.9889	15.4345	0.257119	15.3087	10.6647	2.98963	100.417	7.56672], ...
            [28.1209	1.78503	17.6817	0.15997	13.3335	8.18304	5.14657	20.39	3.70983], ...
            [27.8917	1.73272	18.7614	0.13879	12.6072	7.64412	5.47647	16.8153	2.26001], ...
            [28.9476	1.90182	15.2208	9.98519	15.1	0.261033	3.71601	84.3298	7.97628], ...
            [28.4628	1.68216	18.121	0.142292	12.8429	7.33727	5.59415	16.3535	2.97573], ...
            [29.144	1.83262	15.1726	9.5999	14.7586	0.275116	4.30013	72.029	8.58154], ...
            [28.8131	1.59136	18.4601	0.128903	12.7285	6.76232	5.59927	14.0366	2.39699], ...
            [29.2024	1.77333	15.2293	9.37046	14.5135	0.295977	4.76492	63.3644	9.24354], ...
            [29.1587	1.50711	18.8407	0.116741	12.8268	6.31524	5.38695	12.4244	1.78555], ...
            [29.0818	1.72029	15.43	9.2259	14.4327	0.321703	5.11982	57.056	9.8875], ...
            [29.4936	1.42755	19.3763	0.104621	13.0544	5.93667	5.06412	11.1972	1.01074], ...
            [28.7621	1.67191	15.7189	9.09227	14.5564	0.3505	5.44174	52.0861	10.472], ...
            [28.1894	1.62903	16.155	8.97948	14.9305	0.382661	5.67589	48.1647	11.0005], ...
            [30.419	1.37113	15.2637	6.84706	14.7458	0.165191	5.06795	18.003	6.49804], ...
            [27.3049	1.59279	16.7296	8.86553	15.6115	0.417916	5.83377	45.0011	11.4722], ...
            [30.4156	1.34323	15.862	7.10909	13.6145	0.204633	5.82008	20.3254	8.27903], ...
            [30.7058	1.30923	15.5512	6.71983	14.2326	0.167252	5.53672	17.4911	6.96824], ...
            [27.0059	1.51293	17.7639	8.81174	15.7131	0.424593	5.7837	38.6103	11.6883], ...
            [29.8429	1.32927	16.7224	7.38979	13.2153	0.263297	6.35234	22.9426	9.85329], ...
            [30.9612	1.24813	15.9829	6.60834	13.7348	0.16864	5.92034	16.9392	7.39534], ...
            [16.8819	0.4611	18.5913	8.6216	25.5582	1.4826	5.86	36.3956	12.0658], ...
            [28.0109	1.35321	17.8204	7.7395	14.3359	0.356752	6.58077	26.4043	11.2299], ...
            [30.6886	1.2199	16.9029	6.82872	12.7801	0.212867	6.52354	18.659	9.0968], ...
            [20.6809	0.545	19.0417	8.4484	21.6575	1.5729	5.9676	38.3246	12.6089], ...
            [25.0853	1.39507	18.4973	7.65105	16.8883	0.443378	6.48216	28.2262	12.0205], ...
            [29.5641	1.21152	18.06	7.05639	12.8374	0.284738	6.89912	20.7482	10.6268], ...
            [27.5446	0.65515	19.1584	8.70751	15.538	1.96347	5.52593	45.8149	13.1746], ...
            [21.3985	1.4711	20.4723	0.517394	18.7478	7.43463	6.82847	28.8482	12.5258], ...
            [30.8695	1.1008	18.3481	6.53852	11.9328	0.219074	7.00574	17.2114	9.8027], ...
            [31.0617	0.6902	13.0637	2.3576	18.442	8.618	5.9696	47.2579	13.4118], ...
            [21.7886	1.3366	19.5682	0.488383	19.1406	6.7727	7.01107	23.8132	12.4734], ...
            [32.1244	1.00566	18.8003	6.10926	12.0175	0.147041	6.96886	14.714	8.08428], ...
            [33.3689	0.704	12.951	2.9238	16.5877	8.7937	6.4692	48.0093	13.5782], ...
            [21.8053	1.2356	19.5026	6.24149	19.1053	0.469999	7.10295	20.3185	12.4711], ...
            [33.5364	0.91654	25.0946	0.39042	19.2497	5.71414	6.91555	12.8285	-6.7994], ...
            [34.6726	0.700999	15.4733	3.55078	13.1138	9.55642	7.02588	47.0045	13.677], ...
            [35.3163	0.68587	19.0211	3.97458	9.49887	11.3824	7.42518	45.4715	13.7108], ...
            [35.5631	0.6631	21.2816	4.0691	8.0037	14.0422	7.4433	44.2473	13.6905], ...
            [35.9299	0.646453	23.0547	4.17619	12.1439	23.1052	2.11253	150.645	13.7247], ...
            [35.763	0.616341	22.9064	3.87135	12.4739	19.9887	3.21097	142.325	13.6211], ...
            [35.215	0.604909	21.67	3.5767	7.91342	12.601	7.65078	29.8436	13.5431], ...
            [35.6597	0.589092	23.1032	3.65155	12.5977	18.599	4.08655	117.02	13.5266], ...
            [35.1736	0.579689	22.1112	3.41437	8.19216	12.9187	7.05545	25.9443	13.4637], ...
            [35.5645	0.563359	23.4219	3.46204	12.7473	17.8309	4.80703	99.1722	13.4314], ...
            [35.1007	0.555054	22.4418	3.24498	9.78554	13.4661	5.29444	23.9533	13.376], ...
            [35.8847	0.547751	23.2948	3.41519	14.1891	16.9235	4.17287	105.251	13.4287], ...
            [36.0228	0.5293	23.4128	3.3253	14.9491	16.0927	4.188	100.613	13.3966], ...
            [35.5747	0.52048	22.5259	3.12293	12.2165	12.7148	5.37073	26.3394	13.3092], ...
            [35.3715	0.516598	22.5326	3.05053	12.0291	12.5723	4.7984	23.4582	13.2671], ...
            [34.8509	0.507079	22.7584	2.8903	14.0099	13.1767	1.21457	25.2017	13.1665], ...
            [36.1874	0.511929	23.5964	3.25396	15.6402	15.3622	4.1855	97.4908	13.3573], ...
            [35.7074	0.502322	22.613	3.03807	12.9898	12.1449	5.43227	25.4928	13.2544], ...
            [35.5103	0.498626	22.5787	2.96627	12.7766	11.9484	4.92159	22.7502	13.2116], ...
            [35.0136	0.48981	22.7286	2.81099	14.3884	12.33	1.75669	22.6581	13.113], ...
            [36.5254	0.499384	23.8083	3.26371	16.7707	14.9455	3.47947	105.98	13.3812], ...
            [35.84	0.484938	22.7169	2.96118	13.5807	11.5331	5.66016	24.3992	13.1991], ...
            [35.6493	0.481422	22.646	2.8902	13.3595	11.316	5.18831	21.8301	13.1555], ...
            [35.1736	0.473204	22.7181	2.73848	14.7635	11.553	2.28678	20.9303	13.0582], ...
            [36.6706	0.483629	24.0992	3.20647	17.3415	14.3136	3.49331	102.273	13.3592], ...
            [36.6488	0.465154	24.4096	3.08997	17.399	13.4346	4.21665	88.4834	13.2887], ...
            [36.7881	0.451018	24.7736	3.04619	17.8919	12.8946	4.23284	86.003	13.2754], ...
            [36.9185	0.437533	25.1995	3.00775	18.3317	12.4044	4.24391	83.7881	13.2674], ...
            zeros(1,9)};
            % c_values{68} = zeros(1,9); % disable Ge
            % c_values{110} = zeros(1,9); % disable Sb
            % c_values{113} = zeros(1,9); % disable Te
            
            % match element name without trailing number
            match = regexp(element,'[a-zA-Z+]+','match');
            % index of the elements constants (lower for "non-case-sensitivity")
            ind = find(ismember(lower(c_elements),lower(match{1})));
            % set constants for the element
            if isempty(ind)
                fprintf('Unknown element %s, setting f(q) to zero.\n', match{1})
                c = zeros(1,9);
            else
                c = c_values{ind};
            end

            f = eval(fstr);
        end
        
        
        % convert the build structure objects to something more handy
        % populate obj.zpos{i} and obj.elements{i}
        function convert_struct_list(obj)
            
            for i=1:obj.n_struct
                struct_size = length(obj.struct_list{i}.atompos);
                elements = repmat({''},struct_size,1);
                zpositions = zeros(struct_size,1);
                for s=1:struct_size
                    elements{s} = obj.struct_list{i}.atompos{s}.element;
                    zpositions(s) = obj.struct_list{i}.atompos{s}.pos(3);
                end
                obj.zpos{i} = zpositions;
                obj.elements{i} = elements;
            end
            
        end
        
        % convert the build structure objects to something more handy
        % populate obj.zpos{i} and obj.elements{i}
        function convert_stack_list(obj)

            for i=1:obj.n_struct
                obj.zpos{i} = obj.struct_list{i}.zpos';
                obj.elements{i} = obj.struct_list{i}.elements';
            end            
        end
    end

end