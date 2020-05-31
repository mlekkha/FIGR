function [xntg,tt,nucleusNames,geneNames] = readGeneExprFiles (app)

disp (sigmoid(1.0));

fnameTN = "tn.txt";
fnameXNTG = "xntg.txt";

%======== READ THREE-COLUMN FILE CONTAINING SAMPLE/TIME/NUCLEUS INFO
fid = fopen (fnameTN, 'r');        %Opens File containing time points
matrix = textscan (fid, '%s%s%s', 'Delimiter', '\t', 'HeaderLines', 1, 'MultipleDelimsAsOne', 1);
fclose (fid);
ss = strtrim (matrix{1,1});      % list of sample IDs for each sample s (whitespace removed!)
ts = str2double (matrix{1,2});   % list of times t_s for each sample s
ns = str2double (matrix{1,3});   % list of nuclei n_s for each sample s
tt = unique (ts);  % actual time T_t corresponding to timepoint index t
nn = unique (ns);  % nucleus label N_n corresponding to nucleus index n
smax = numel (ss); % number of samples
tmax = numel (tt); % number of unique timepoints
nmax = numel (nn); % number of unique nuclei

%======== COMPARE SAMPLE NAMES BETWEEN TWO FILES
fid = fopen(fnameXNTG , 'r');
dummy = fgetl (fid);
dummy = split (strtrim (dummy)); % split into strings
ssFromMainFile = dummy(2:end);
assert (numel(ssFromMainFile) == smax, 'ERROR: incompatible numbers of samples in data files!');
assert (isequal (ssFromMainFile, ss), 'ERROR: differing sample names in data files!');
%======== READ EACH LINE OF THE FILE AND PARSE IT INTO THE ARRAY
geneNames = [];                          % could be "gg"
x_gs = [];
while true
dummy = fgetl (fid);                 % read one line from the file (as char vec)
if ~ischar(dummy); break; end
dummy = split (strtrim (dummy));     % split line into strings
geneNames = [geneNames ; dummy(1)];  % first item is gene name
dummy = str2double ( dummy(2:end) )'; % remaining items are gene expression levels
x_gs = [x_gs ; dummy];
end
fclose (fid);
geneNames = string (geneNames);          % convert from CELL ARRAY OF CHAR VECTORS to ARRAY OF STRINGS
gmax = numel (geneNames);

%======== AVERAGE OVER REPLICATES
[tnu, ia, ic] = unique ([ts ns], 'rows', 'stable');
umax = numel (ia);  % number of unique datapoints
x_gu = zeros ([gmax umax]);
for g = 1:gmax
x_gu(g,:) = accumarray (ic, x_gs(g,:), [], @mean)';
end

%======== RESHAPE x(n,t,g) INTO x(g,u)
xntg = zeros ([nmax tmax gmax]);
for u = 1:umax
tn = tnu(u,:);  % time and position corresponding to sample s
t = tn(1);
n = tn(2);
indext = find (tt==t);  % time index corresponding to sample s
indexn = find (nn==n);
xntg (indexn,indext,:) = x_gu (:,u);
end

%======== FINALIZE RETURN-VALUES
nucleusNames = string (nn);
end
