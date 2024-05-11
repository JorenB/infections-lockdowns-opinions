function cfg = cfgRead(filename)
	jscfg = fileread(filename);

	cfg = jsondecode(jscfg);
end
