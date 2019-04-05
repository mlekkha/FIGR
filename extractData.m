function [numGenes,numExternals,grn,Xntg,nuclei,tt] = ...
    extractData (opts, inFileName)

debug = opts.debug;

% READS:   inFileName (input file name, e.g., 'hkgn58c14k1_002')
% RETURNS: Rg      (1:numGenes)
%          lambdag (1:numGenes)
%          Tgg     (1:numGenes, 1:numGenes+numExternals+1)
%          Xntg    (1:numNuclei, 1:numTimepoints, 1:numGenes+numExternals)
%          nuclei  (1:numNuclei)
%          tt      (1:numTimepoints)

fprintf ('======== extractData() ===============================\n');

%======== PREPARE INPUTS INCLUDING EXTERNAL INPUTS ========
fid = fopen (inFileName, 'r');

findTag (fid, '$problem');
numGenes = readTaggedMatrix (fid, 'number_of_genes');
numExternals = readTaggedMatrix (fid, 'number_of_external_inputs') + 1;

findTag (fid, '$eqparms');
Rg = readTaggedMatrix (fid, 'promoter_strengths'); Rg=Rg(:); % make into column vector1
Tgg = readTaggedMatrix (fid, 'genetic_interconnect_matrix');
Egg = readTaggedMatrix (fid, 'external_input_strengths');
Mg = readTaggedMatrix (fid, 'maternal_connection_strengths'); Mg=Mg(:);
hg = readTaggedMatrix (fid, 'promoter_thresholds');           hg=hg(:);
Dg = readTaggedMatrix (fid, 'diffusion_parameter(s)');        Dg=Dg(:);
halfLives = readTaggedMatrix (fid, 'protein_half_lives');
lambdag = log(2)./halfLives;                        lambdag=lambdag(:);

%-------- READ Hunchback, Kruppel, Giant, Knirps
% Lineage, Time, Hb, Kr, Gt, Kni
findTag (fid, '$bias_wt'); % (58 x 1) x 6
LTHKGN0 = readTaggedMatrix (fid, 'lin');  % pretend that the line containing 'lin' was the tag
findTag (fid, '$facts_wt'); % (58 x 8) x 6
LTHKGNt = readTaggedMatrix (fid, 'lin');
LTHKGN = [LTHKGN0 ; LTHKGNt];  % (58 x 9) x 6
%-------- READ Bicoid
LB = readTaggedMatrix (fid, '$bcd_wt');  % 58 x 2
%-------- READ Caudal, Tailless
findTag (fid, '$ext_wt');      % (58 x 9) x 6
LTCT = readTaggedMatrix (fid, 'lin');
LTCT = LTCT(  LTCT(:,1)>8191,  :  );    % SELECT ONLY CYCLE-14 DATA

%-------- REARRANGE INTO Xntg(nucleus,timepoint,gene) FORMAT
nuclei = unique (LTCT(:,1)); % [8227 8228 ... 8249]'
tt     = unique (LTCT(:,2)); % [0 3.125 ... 46.875]'
%nuclei        = unique (LTHKGN(:,1)); % [8227 8228 ... 8249]'
%tt            = unique (LTHKGN(:,2)); % [0 3.125 ... 46.875]'
%nuclei = LB(:,1);
numNuclei     = numel (nuclei);              % 58
numTimepoints = numel (tt);                  % 9

HKGN = LTHKGN (:,3:end);              % (58 x 9) x 4
B    = LB(:,2);                       % 58 x 1
B    = repmat (B, numTimepoints, 1);  % (58 x 9) x 1
CT   = LTCT (:,3:end);                % (58 x 9) x 2
Xntg = [HKGN  B  CT];                 % (58 x 9) x (4+1+2)
Xntg = reshape (Xntg, numNuclei, numTimepoints, []);


%-------- Tgg(gene,gene)
% Columns are: Hb, Kr, Gt, Kni, Bcd, Cad, Tll, and lastly h_i (threshold)
Tgg = [Tgg  Mg  Egg];

grn = {};
grn.Tgg = Tgg;
grn.hg = hg;
grn.Rg = Rg;
grn.lambdag = lambdag;
grn.Dg = Dg;

if debug
	fprintf ('======== ARRANGING DATA ==============\n');
	fprintf ('numNuclei          = %d \n',numNuclei);
	fprintf ('numTimepoints      = %d \n',numTimepoints);
	fprintf ('numGenes           = %d \n',numGenes);
	fprintf ('numExternals       = %d \n',numExternals);
	fprintf ('size(HKGN) = %d x %d\n', size(HKGN,1), size(HKGN,2));
	fprintf ('size(B)    = %d x %d\n', size(B,1), size(B,2));
	fprintf ('size(CT)   = %d x %d\n', size(CT,1), size(CT,2));
	fprintf ('size(Xntg) = %d x %d x %d \n', size(Xntg));
	fprintf ('tt     = \n %s\n'    , mat2str_compact (tt') );
	fprintf ('nuclei = \n %s\n\n', num2str (nuclei') );
	fprintf ('Rg = \n'); disp (Rg);
	fprintf ('Tgg = \n'); disp (Tgg);
	fprintf ('Dg = \n'); disp (Dg);
	fprintf ('half-lives = \n'); disp (halfLives);
	fprintf ('lambdag = \n'); disp (lambdag);
end

fclose (fid);






	%======== NESTED FUNCTIONS FOR STRING I/O ===============
    function out = readTaggedMatrix (fid, tag)
        findTag (fid, tag);
        out = readNextMatrix (fid);
        if debug
       		if numel(out) < 20
        	    disp (out);
        	else
            	fprintf ('      [%d x %d array]\n', size(out,1), size(out,2));
        	end
        end
    end

	% Read lines one by one from file "fid" until a line containing string "tag" is found
    function [] = findTag (fid, tag)  %WARNING: THIS ONLY SEARCHES IN FORWARD DIRECTION!
        while ~feof(fid)
            tline = fgetl(fid);
            if (~isempty(strfind (tline, tag)))
                if debug; fprintf ('Tag: %s \n', tag); end;
                break;
            end
        end
    end
    function out = readNextVector (fid)
        out = fscanf (fid, '%f', inf); % Return a col vec regardless of arrangement of nums in file
    end
    function out = readNextMatrix (fid)
        while true
            %==== Remember current position
            currentFilePos = ftell (fid);
            %==== Count number of numbers on next line
            tline = fgets(fid);  % include newline so textscan doesn't choke on empty string
            data = textscan (tline, '%f');
            nColumns = numel(data{1});
            if nColumns > 0; break; end
        end
        %==== Rewind to where we were
        fseek (fid, currentFilePos, 'bof');
        %==== Now read data for real
        %data = textscan (fid, '%f', inf, 'Delimiter', {' ','\b','\t'});
        data = textscan (fid, '%f', inf);
        data = data{1};
        out = reshape (data, nColumns, [])';
    end
    function out = mat2str (mat)     % Convert matrix to string
        fmt = [ repmat('%g\t', 1, size(mat, 2))     '\n' ];
        out = sprintf (fmt, transpose(mat));
    end
    function out = mat2str_compact (mat)
        fmt = [ repmat('%g ', 1, size(mat, 2))     '\n' ];
        out = sprintf (fmt, transpose(mat));
    end
    function out = multiline (strs)  % Construct multiline string literal
        strs(:,2) = {'\n'};
        strs=strs';
        out = cat(2, strs{:});
    end
    function out = insertBlanksInString (str)
        out = sprintf ('%8c', str);
    end
    function out = flatten (array)
        out = reshape (array, 1, numel(array));
    end
    %======= END NESTED FUNCTIONS
end