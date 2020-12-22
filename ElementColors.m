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

function [ret] = ElementColors(name, num)

switch lower(name)
   % Henning/Marc Diss Colors
   case {'sb','antimony'}
      palette = {[255 59 101] [230 0 50] [184 0 40]};
   case {'ge','germanium'}
      palette = {[98 145 205] [52 101 164] [41 80 131]}; 
   case {'ti','titanium'}
      palette = {[156 236 78] [115 210 22] [92 168 17]}; % green
   case {'te','tellurium'}
      palette = {[223 223 223] [192 192 192] [153 153 153]};
   case {'sl','csl'}
      palette = {[185 0 185] [128 0 128] [100 0 100]}; % purple
   case {'ga'}
      palette = {[16 51 153] [24 76 230] [101 136 238]}; % blue
   case {'al'}
      palette = {[16 51 153] [24 76 230] [101 136 238]}; % blue
   case {'as'}
      palette = {[254 251 170] [253 248 86] [252 244 2]}; % yellow
   case {'in'}
      palette = {[253 173 110] [252 113 2] [216 96 2]}; % orange
%    case {'te'}
%       palette = {[252 233 79] [237 212 0] [196 160 0]};
%    case {'sb'}
%       palette = {[252 175 62] [245 121 0] [206 92 0]};
   case {'brown'}
      palette = {[233 185 110] [193 125 17] [143 89 2]};
   case {'green'}
      palette = {[138 226 52] [115 210 22] [78 154 6]};
   case {'blue'}
      palette = {[114 159 207] [52 101 164] [32 74 135]}; 
%    case {'ti'}
%       palette = {[114 159 207] [52 101 164] [32 74 135]}; 
%    case {'ge'}
%       palette = {[173 127 168] [117 80 123] [92 53 102]}; 
   case {'red'}
      palette = {[239 40 40] [204 0 0] [164 0 0]}; 
   case {'aluminum'}
      palette = {[211 215 207]};   
   case {'darkaluminum'}
      palette = {[85 87 83]};       
   case {'sky blue'}
      palette = {[52 101 164]};         
   case {'gray'}
      palette = {[238 238 236] [211 215 207] [186 189 182] ...
                 [136 138 133] [85 87 83] [46 52 54] [0 0 0]};
    case {'vdw'}
      palette = {[238 238 236] [211 215 207] [186 189 182] ...
                 [136 138 133] [85 87 83] [46 52 54] [0 0 0]};
    case {'vacancy'}
      palette = {[238 238 236] [211 215 207] [186 189 182] ...
                 [136 138 133] [85 87 83] [46 52 54] [0 0 0]};
    case {'gap'}
      palette = {[238 238 236] [211 215 207] [186 189 182] ...
                 [136 138 133] [85 87 83] [46 52 54] [0 0 0]};
    case {'black'}
        palette = {[0 0 0]};
    otherwise
      palette = {[238 238 236] [211 215 207] [186 189 182] ...
          [136 138 133] [85 87 83] [46 52 54]};
      disp('Unknown Colorpalette, using "gray"');
end
   nops = size(palette);
   if num > nops(2)
      ret = [0 0 0];
      fprintf('Unknown Color ElementColors(%s,%i)',name, num);      
   else   
     ret = palette{num};
   end  
   ret = ret/255;
end