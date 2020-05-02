classdef noname < SerialLink
 
	properties
	end
 
	methods
		function ro = noname()
			objdir = which('noname');
			idx = find(objdir == filesep,2,'last');
			objdir = objdir(1:idx(1));
			 
			tmp = load(fullfile(objdir,'@noname','matnoname.mat'));
			 
			ro = ro@SerialLink(tmp.sr);
			 
			 
		end
	end
	 
end
