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
 %}

classdef buildStack < handle
    properties
       zpos = [];
       space = 0;
       next_zpos = 0;
       elements = {};
       z = 0;
    end
    
    methods

        % place to initialize things
        function obj = buildStack(varargin)
            %
        end
        
        function setSpace(obj)
            obj.space = obj.zpos(end);
            obj.next_zpos = obj.space;
            obj.z = obj.zpos(end);
        end
            
        function addAtom(obj,elements,zpos)
            if isa(elements,'cell')
                obj.elements = [obj.elements elements];
            elseif isa(elements,'char')
                obj.elements{end+1} = elements;
            end
            obj.zpos = [obj.zpos obj.space+zpos];
            obj.setSpace
        end
        
        function multiaddAtom(obj,elements,zpos,space,n)
            for i=1:n
                obj.addAtom(elements,zpos);
                obj.addDistance(space);
            end
        end
        
        function multiaddBlock(obj,block,n)
            for i=1:n
                obj.addBlock(block);
            end
        end
        
        function addBlock(obj,block)
            obj.addAtom(block.elements,block.zpos);
            obj.addDistance(block.space);
        end
        
        function addDistance(obj,distance)
            obj.space = obj.space + distance;
            obj.z = obj.z + distance;
        end
        
        function closeStack(obj)
            obj.z = obj.space;
        end
        
        % https://materials.springer.com/isp/crystallographic/docs/sd_0260901
        function block = SnTeBL111(~,c) % SnTe(111) blocks
            % a_SnTe = 6.301;
            % c_SnTe = sqrt(3)*a_SnTe = 10.9137
            % a_SnTe/sqrt(3) = SnTe(111) blocks
            block.elements = {'Te', 'Sn'};
            block.zpos = [0.0000 1/6].*c;
            block.space = 1/6*c;
        end     
        
        % https://materials.springer.com/isp/crystallographic/docs/sd_0260901
        function block = SnTeBL100(~,c) % SnTe(111) blocks
            % a_SnTe = 6.301;
            % c_SnTe = sqrt(3)*a_SnTe = 10.9137
            % a_SnTe/sqrt(3) = SnTe(111) blocks
            block.elements = {'Te', 'Sn'};
            block.zpos = [0.0000 0.500].*c;
            block.space = 0.500*c;
        end
        
        function block = GeTeBL111(~,c)   
            block.elements = {'Te', 'Ge'};
            block.zpos = [0.0000 0.1419].*c;
            block.space = 0.1915*c;    
        end
        
        function block = TGG(~,c)   
            block.elements = {'Te', 'Ge', 'Ge'};
            block.zpos = [0.0000 0.1419 0.2869].*c;
            block.space = 0.1419*c;    
        end
        
        function block = vacuumGeTeBL111(~,c)   
            block.elements = {'vacuum', 'vacuum'};
            block.zpos = [0.0000 0.1419].*c;
            block.space = 0.1915*c;    
        end
        
        % alias for GeTeBL111
        function block = BL(obj,c)
            block = obj.GeTeBL111(c);
        end
        
        % alias for GeTeBL111
        function block = GeTeBL(obj,c)
            block = obj.GeTeBL111(c);
        end
        
        % alias for SnTeBL111
        function block = SnTeBL(obj,c)
            block = obj.SnTeBL111(c);
        end
        
        function block = TeTeTi(~,c)
            %c_TiTe2 = 6.49800;
            block.elements = {'Te', 'Te', 'Ti'};
            block.zpos = [0    0.4744    0.7372].*c;
            block.space = 0.26280*c;
        end
        
        function block = GaAs(~,c)
            %c_GaAs = 5.65370;
            block.elements = {'Ga', 'As'};
            block.zpos = [0    0.25].*c;
            block.space = 0.25*c;
        end
        
        function block = AlAs(~,c)
            %c_AlAs = 5.66080;
            block.elements = {'Al', 'As'};
            block.zpos = [0    0.25].*c;
            block.space = 0.25*c;
        end
        
        function block = InAs(~,c)
            %c_InAs = 6.04000;
            block.elements = {'In', 'As'};
            block.zpos = [0    0.25].*c;
            block.space = 0.25*c;
        end
        
        function block = TiTe2(~,c)
            %c_TiTe2 = 6.49800;
            block.elements = {'Ti', 'Te', 'Te'};
            block.zpos = [0    0.26280    0.7372].*c;
            block.space = 0.2628*c;
        end
        
        % quintuple structure with inserted Sb2 forming a Sb8Te9
        % DOI: 10.1039/b500695c
        function block = QT_Sb2(~,c)
            % c_Sb2Te3 = 30.6069; c_Sb8Te9 = 102.6900;
            %[ 0    2.0537    3.6881 6.1773 7.7248 10.2150 11.8494]
            block.elements = {'Te', 'Sb', 'Te', 'Sb', 'Sb', 'Te', 'Sb'};
            block.zpos = [0    0.0671    0.1205    0.2018    0.2524    0.3337 0.3871]*c;
            block.space = 0.0671*c;
        end
        
        % quintuple structure to build a Momand/Kooi CSL
        function block = GST124(~,c)
            % c_GST124 = 41.68600;
            block.elements = {'Te', 'Sb', 'Te', 'Te', 'Sb', 'Te', 'Ge'};
            block.zpos = [0    0.0504    0.0895    0.1572    0.1964    0.2468    0.2901]*c;
            block.space = 0.0433*c;
        end
        
        % quintuple structure to build a Momand/Kooi CSL
        function block = QT_Kooi(~,c)
            % c_Sb2Te3 = 30.6069;
            block.elements = {'Te', 'Sb', 'Te', 'Te', 'Sb'};
            block.zpos = [0.0000 0.0671 0.1205 0.2128 0.2662]*c;
            block.space = 0.0671*c;
        end
        
        % quintuple structure to build a Momand/Kooi CSL
        function block = vacuumQT_Kooi(~,c)
            % c_Sb2Te3 = 30.6069;
            block.elements = {'vacuum', 'vacuum', 'vacuum', 'vacuum', 'vacuum'};
            block.zpos = [0.0000 0.0671 0.1205 0.2128 0.2662]*c;
            block.space = 0.0671*c;    
        end
        
        function block = QL_Kooi(obj,c)
            block = obj.QT_Kooi(c);
        end
        
        % quintuple block with variable vdW distance
        function block = QT_Kooi_vdw(~,c,d_vdw)
            % d_vdw=2.8250;
            block.elements = {'Te', 'Sb', 'Te', 'Te', 'Sb'};
%             block.zpos = [[0.0000 0.0671 0.1205].*c [0.2128 0.2662].*c+d_vdw];
            block.zpos = [[0.0000 0.0671 0.1205].*c [0.1205 0.1739].*c+d_vdw];
            block.space = 0.0671*c;
        end
        
        % quintuple block with variable vdW distance, keep c constant
        function block = QT_Kooi_vdw_c_const(~,c,d_vdw)
            l = 0.0671; % long bond
            s = 0.0534; % shot bond
            ratio = 1/(l+s);
            l_s=0.5*((1/3-d_vdw/c)); % bond length of l+s
            l = l_s*l*ratio; % stretched long bond
            s = l_s*s*ratio; % stretched short bond
            
            % check (yields 1/3)
            2*(l+s)+d_vdw/c;
            
            block.elements = {'Te', 'Sb', 'Te', 'Te', 'Sb'};
            block.zpos = [0 0 0 0 0];
            bonds = [l s d_vdw/c s l].*c;
            for i=2:length(block.elements)
               block.zpos(i) = (block.zpos(i-1)+bonds(i-1));
            end
            block.space = bonds(end);
        end
        
        
        % QT block with vdW gaps at the end
        function block = QT(~,c)
            block.elements = {'Te', 'Sb', 'Te', 'Sb', 'Te'};
            block.zpos = [0.0000 0.0534 0.1205 0.1876 0.2410]*c;
            % vdW gap at the end
            block.space = 0.0923*c;
        end
        
        % QL block with vdW gaps at the end
        function block = QL(obj,c)
            block = obj.QT(c);
        end
        
        % QT block with vdW gaps at the end
        function block = Sb2Te3(~,c)
            block.elements = {'Te', 'Sb', 'Te', 'Te', 'Sb'};
            block.zpos = [0.0000 0.06547 0.12053 0.21280 0.26787]*c;
            % vdW gap at the end
            block.space = (0.33333-0.26787)*c;
        end

    end
    
    methods(Static)
 
    end
end
