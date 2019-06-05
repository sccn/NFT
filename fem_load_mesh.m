% Usage: [Coord, Elem, Sig] = load_mesh (FILENAME)

% Zeynep Akalin Acar, 2015
%
% Load a FEM mesh file from FILENAME.
% supports old-style and parallel mesh files
% with linear (8 node) and quadratic (20 node) elements.

function [Coord, Elem, Sig] = fem_load_mesh (fname)
	fid = fopen(fname, 'r');
	if fid == -1
		error ('Failed to open %s\n', fname);
	end

	line = fgets(fid);
	[nn, cnt] = sscanf(line, '%d');
	if (cnt == 0 || nn < 0)
		error('Invalid mesh file\n');
	end

	if (nn == 0)
		line = fgets(fid);
		[nn, cnt] = sscanf(line, '%d', 1);
		if cnt == 0 || nn ~= 0
			error ('Invalid mesh file\n');
		end
		line = fgets(fid);
		[nn, cnt] = sscanf(line, '%d', 1);
		if cnt == 0 || nn <= 0
			error ('Invalid mesh file\n');
		end
	end

	[C, cnt] = fscanf(fid, '%d %g %g %g', [4, nn]);
	if cnt ~= nn*4
		error ('Error reading nodes\n');
	end

	idx = 1:nn;
	if idx ~= C(1,:)
		error ('Error in node indices\n');
	end

	line = fgets(fid);	% skip end of line
	line = fgets(fid);
	[tmp, cnt] = sscanf(line, '%d');
	if cnt <= 0
		error('Invalid mesh file\n');
	end

	ns = 20;
	ne = tmp(1);
	if (cnt == 2)
		ns = tmp(2);
	end
	
	if (ns ~= 20 && ns ~= 8 && ns~= 4 && ns ~= 10)
		error('Invalid nodes/element');
	else
		[E, cnt] = fscanf(fid, '%d', [ns + 1, ne]);
	end

	if cnt ~= ne * (ns+1)
		error('Error reading elements\n');
	end

	idx = 1:ne;
	if idx ~= E(1,:)
		error('Error in element indices\n');
	end

	line = fgets(fid);	% skip end of line

	line = fgets(fid);
	[tmp, cnt] = sscanf(line, '%d');
	S = zeros(ne, 2);

	if (cnt == 0)
		warning('No conductivity information\n');
	elseif (tmp == nn)
		warning('Node conductivity not supported\n');
	elseif (tmp == ne)
		pos = ftell(fid);
		[S, cnt] = fscanf(fid, '%d %g', [2, ne]);
		if cnt ~= ne*2
			fseek(fid, pos, 'bof');
			[S, cnt] = fscanf(fid, '%d C%d', [2, ne]);
			if cnt ~= ne*2
				error('Failed to read sigma\n');
			end
		end
		if idx ~= S(1,:)
			error('Error in conductivity indices\n');
		end
		line = fgets(fid);	% skip end of line
	end

	Coord = C(2:4,:)';
	Elem  = E(2:ns+1,:)';
	Sig   = S(2,:)';

	fclose(fid);

end
