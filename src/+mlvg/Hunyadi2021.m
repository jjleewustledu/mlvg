classdef Hunyadi2021 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
	%% HUNYADI2021 is derived from Levente Hunyadi's methods for B-splines from
    %  https://www.mathworks.com/matlabcentral/fileexchange/27374-b-splines?s_tid=srchtitle
    %  using licensing detailed in license_Hunyadi.txt.

	%  $Revision$
 	%  was created 25-Mar-2021 21:25:12 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlvg/src/+mlvg.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = Hunyadi2021(varargin)
            % Illustrates B-spline curve estimation without knowing parameter values.
            % Copyright 2010 Levente Hunyadi            
            
            addpath(genpath(fullfile(getenv('HOME'), 'MATLAB-Drive', 'bspline')))
            
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

