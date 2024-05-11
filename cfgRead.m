% Copyright 2024 Alexandra Teslya & Joren Brunekreef
function cfg = cfgRead(filename)
	jscfg = fileread(filename);

	cfg = jsondecode(jscfg);
end
