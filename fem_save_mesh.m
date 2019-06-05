% Usage: save_mesh (FILENAME, Coord, Elem, Sig)
%
% zeynep Akalin Acar 2015
% Saves a FEM mesh file to FILENAME.
% saves old-style mesh files
% supports linear (8 node) and quadratic (20 node) hexahedral elements.
% and linear (4 node) and quadratic (10 node) tetrahedral elements.

function fem_save_mesh (fname, Coord, Elem, Sig)
	fid = fopen(fname, 'w');
	if fid == -1
		error ('Failed to open %s\n', fname);
	end
	nn = size(Coord, 1);
	if size(Coord, 2) ~= 3
		error('invalid coordinate matrix');
	end
	ne = size(Elem, 1);
	nne = size(Elem, 2);

	if (nne ~= 8 && nne ~= 20 && nne ~= 4 && nne ~= 10)
		error('Invalid node/elem');
	end

	ns = size(Sig, 1);
	if ns ~= ne
		error('Invalid conductivity size');
	end
	if size(Sig, 2) ~= 1
		error('invalid conductivity vector');
	end

	fprintf(fid, '%d\n', nn);

	C = zeros(4, nn);
	C(1,:) = 1:nn;
	C(2:4,:) = Coord';
	fprintf(fid, '%d %g %g %g\n', C);

	fprintf(fid, '%d %d\n', ne, nne);

	if (nne == 4)
		E = zeros(5, ne);
		E(1,:) = 1:ne;
		E(2:5,:) = Elem';
		fprintf(fid, '%d %g %g %g %g\n', E);
	elseif (nne == 8)
		E = zeros(9, ne);
		E(1,:) = 1:ne;
		E(2:9,:) = Elem';
		fprintf(fid, '%d %g %g %g %g %g %g %g %g\n', E);
	elseif (nne == 10)
		E = zeros(11, ne);
		E(1,:) = 1:ne;
		E(2:11,:) = Elem';
		fprintf(fid, '%d %g %g %g %g %g %g %g %g %g %g\n', E);
    else
		E = zeros(21, ne);
		E(1,:) = 1:ne;
		E(2:21,:) = Elem';
		fprintf(fid, '%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n', E);
    end

	fprintf(fid, '%d\n', ns);
    if (min(Sig) == 1) && (max(Sig - floor(Sig))==0)
        % Sig contains classes
        for i = 1:ns
            fprintf(fid, '%d C%d\n', i, Sig(i));
        end
    else
        % Sig contains conductivity values
    	S = zeros(2, ns);
    	S(1,:) = 1:ns;
    	S(2,:) = Sig';
    	fprintf(fid, '%d %g\n', S);
    end
    
	fclose(fid);
end
