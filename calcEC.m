%Calculate_evolutionary_constraints.m  (Version 2.0)
%
%This software serves to extract evolutionary couplings from protein multiple sequence alignments.
%
%Contributors
%
%Version 1.0  (2011) - Andrea Pagnani, Martin Weigt
%Version 2.0  (2011) - Lucy Colwell, Debora Marks, Rob Sheridan
%
%Memorial Sloan-Kettering Cancer Center - Group of Chris Sander
%Harvard Medical School - Group of Debora Marks
%Politecnico di Torino - Group of Riccardo Zecchina
%
%References
%
%Marks et al - 2011
%Morcos et al - 2011
%Hopf et al - 2012
%
%License and distribution
%
%Free for academic, non-commercial use.  No redistribution.
%Distributed to registered users. Registration for update purposes.
%Registration on www.evfold.org - or send email to Rob Sheridan (sheridan@cbio.mskcc.org)
%Cite the above references in publications using this code.

function calcEC(msa_fasta_filename, seqid_of_interest, outputfile)
% usage: calculate_evolutionary_constraints('PF00071_v25.fa', 'RASH_HUMAN', 'PF00071_P01112_MI_DI.txt')
theta = getenv('DI_THETA');
if (size(theta,2) == 0)
    theta = 0.3;
else
    theta = str2num(theta);
end
pseudocount_weight = getenv('DI_PSEUDOCOUNT_WEIGHT');
if (size(pseudocount_weight,2) == 0)
    pseudocount_weight = 0.5;
else
    pseudocount_weight = str2num(pseudocount_weight);
end
[Pij_true, Pi_true, alignment_width, q, encoded_seq_of_interest, focus_to_uniprot_offset_map] = read_alignment(msa_fasta_filename, seqid_of_interest, theta);
% Pij_true is the weighted MI.
[Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q);
number2letter_map = create_number2letter_map();
C = Compute_C(Pij, Pi, alignment_width, q);
invC = inv(C);
clear C;
fp = fopen(outputfile, 'w');
for i=1:(alignment_width-1)
    for j=(i+1):alignment_width
        % mutual information
        [MI_true, ~, ~] = calculate_mi(i, j, Pij_true, Pi_true, q);
        fprintf(fp, '%d %s %d %s %g ', focus_to_uniprot_offset_map(i), number2letter_map(encoded_seq_of_interest(i)), focus_to_uniprot_offset_map(j), number2letter_map(encoded_seq_of_interest(j)), MI_true);
        % direct information from mean field
        W_mf = ReturnW(invC, i, j, q);
        % note this is of invC, so we are computing
        % exp(-(invC_ij) for each amino acid pair.
        DI_mf_pc = bp_link(i, j, W_mf, Pi, q);
        fprintf(fp, '%g ', DI_mf_pc);
        fprintf(fp, '\n');
    end
end
fclose(fp);
end

function [Pij_true, Pi_true, alignment_width, q, encoded_seq_of_interest, focus_to_uniprot_offset_map] = read_alignment(msa_fasta_filename, seqid_of_interest, theta)
[encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest);
encoded_seq_of_interest = encoded_focus_alignment(focus_index_of_interest,:);
[alignment_height,alignment_width] = size(encoded_focus_alignment);
W = ones(1, alignment_height);
if(theta > 0.0)
    W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm')<theta))));
end
Meff=sum(W);
q = max(max(encoded_focus_alignment));
Pij_true = zeros(alignment_width, alignment_width, q, q);
Pi_true = zeros(alignment_width, q);
for j=1:alignment_height
    for i=1:alignment_width
        Pi_true(i, encoded_focus_alignment(j, i)) = Pi_true(i, encoded_focus_alignment(j, i)) + W(j);
    end
end
Pi_true = Pi_true/Meff;
for l=1:alignment_height
    for i=1:alignment_width-1
        for j=i+1:alignment_width
            Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) + W(l);
            Pij_true(j, i, encoded_focus_alignment(l, j), encoded_focus_alignment(l, i)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j));
        end
    end
end
Pij_true = Pij_true/Meff;

scra = eye(q, q);
for i=1:alignment_width
    for alpha=1:q
        for beta=1:q
            Pij_true(i, i, alpha, beta) = Pi_true(i, alpha) * scra(alpha, beta);
        end
    end
end
end

function s2=cellstr2cell(s)
s=strtrim(s);
s2={};
if s(1)=='{'
    s=s(2:end-1);
    idx=find(s==',');
    bidx=1;
    for k=1:length(idx)
        s2{k}=s(bidx:idx(k)-1);
        bidx=idx(k)+1;
    end
    s2{end+1}=s(bidx:end);
else
    s2=s;
end
end

function [encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest)
METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = getenv('DI_METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES'); % 1 = change them to gaps .. 2 = mask entire sequence
if (size(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES,2) == 0)
    METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = 2;
else
    METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = str2num(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES)
end
%Changes here
randomize_order=1; %randomize the order (except the first) for both full_alignment and f2
if ischar(msa_fasta_filename)
    msa_fasta_filename = cellstr2cell(msa_fasta_filename);
end
if ~iscell(msa_fasta_filename)
    msa_fasta_filename={msa_fasta_filename};
end
%celldisp(msa_fasta_filename)
AlignIdx=zeros(1,length(msa_fasta_filename));
full_alignment = myfastaread(msa_fasta_filename{1});
AlignIdx(1)=length(full_alignment(1).Sequence);
for idx=2:length(msa_fasta_filename)
    f2=myfastaread(msa_fasta_filename{idx});
    mlen=min(length(full_alignment),length(f2));
    %randomize the order (except the first) for both full_alignment and f2
    if randomize_order
        full_alignment=full_alignment([1 randperm(length(full_alignment)-1)+1]);
        f2=f2([1 randperm(length(f2)-1)+1]);
    end
    disp(['Using only ' num2str(mlen) ' Sequences from the alignmet'])
    full_alignment=full_alignment(1:mlen);
    f2=f2(1:mlen);
    for idx2=1:mlen
        full_alignment(idx2).Sequence=[full_alignment(idx2).Sequence f2(idx2).Sequence];
    end
    AlignIdx(idx)=length(full_alignment(1).Sequence);
end
%Changes end
alignment_width = size(full_alignment(1).Sequence, 2);
alignment_height = size(full_alignment, 1);
letter2number_map = create_letter2number_map();
[full_index_of_interest, range_of_interest_start, range_of_interest_end] = find_seq_of_interest(full_alignment, seqid_of_interest);
encoded_focus_alignment = [];
skipped_sequence_counter = 0;
[focuscolumnlist, focus_to_uniprot_offset_map] = scan_sequence_of_interest_for_focus_columns(full_alignment(full_index_of_interest).Sequence, range_of_interest_start, letter2number_map);
for full_alignment_index=1:alignment_height
    focus_alignment_row = full_alignment(full_alignment_index).Sequence(focuscolumnlist);
    encoded_focus_alignment_row = letter2number_map(focus_alignment_row);
    if (size(find(encoded_focus_alignment_row == 0),2) > 0)
        error(['Error: sequence in alignment has illegal characters: ' full_alignment(full_alignment_index).Sequence]);
    end
    if (size(find(encoded_focus_alignment_row <= -2),2) > 0)
        error(['Error: sequence in alignment has dot or lowercase in conserved position: ' full_alignment(full_alignment_index).Sequence]);
    end
    if (size(find(encoded_focus_alignment_row == -1),2) > 0)
        if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 1)
            encoded_focus_alignment_row(find(encoded_focus_alignment_row == -1)) = 1;
        else
            if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 2)
                continue %skip sequences with ambiguous residues
            else
                error('Internal Error');
            end
        end
    end
    encoded_focus_alignment(size(encoded_focus_alignment,1) + 1,:) = encoded_focus_alignment_row;
    if (full_alignment_index == full_index_of_interest)
        focus_index_of_interest = size(encoded_focus_alignment,1);
    end
end
end

function [Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(alignment_width, alignment_width, q, q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(alignment_width, q);
scra = eye(q);
for i=1:alignment_width
    for alpha = 1:q
        for beta = 1:q
            Pij(i, i, alpha, beta) = (1.-pseudocount_weight)*Pij_true(i, i, alpha, beta) + pseudocount_weight/q*scra(alpha, beta);
        end
    end
end
end

function C = Compute_C(Pij, Pi, alignment_width, q)
C=zeros(alignment_width*(q-1), alignment_width*(q-1));
for i=1:alignment_width
    for j=1:alignment_width
        for alpha=1:q-1
            for beta=1:q-1
                C(mapkey(i, alpha, q), mapkey(j, beta, q)) = Pij(i, j, alpha, beta) - Pi(i, alpha)*Pi(j, beta);
            end
        end
    end
end
end

function A=mapkey(i, alpha, q)
A = (q-1)*(i-1)+alpha;
end

function [M, s1, s2] = calculate_mi(i, j, P2, P1, q)
M = 0.;
for alpha=1:q
    for beta = 1:q
        if(P2(i, j, alpha, beta)>0)
            M = M + P2(i, j, alpha, beta)*log(P2(i, j, alpha, beta) / P1(i, alpha)/P1(j, beta));
        end
    end
end

s1=0.;
s2=0.;
for alpha=1:q
    if(P1(i, alpha)>0)
        s1 = s1 - P1(i, alpha) * log(P1(i, alpha));
    end
    if(P1(j, alpha)>0)
        s2 = s2 - P1(j, alpha) * log(P1(j, alpha));
    end
end
end

function W=ReturnW(C, i, j, q)
W = ones(q, q);
W(1:q-1, 1:q-1) = exp(-C(mapkey(i, 1:q-1, q), mapkey(j, 1:q-1, q)));
end

function DI = bp_link(i, j, W, P1, q)
[mu1, mu2] = compute_mu(i, j, W, P1, q);
DI = compute_di(i, j, W, mu1, mu2, P1);
return;
end

function [mu1, mu2] = compute_mu(i, j, W, P1, q)
epsilon=1e-4;
diff =1.0;
mu1 = ones(1, q)/q;
mu2 = ones(1, q)/q;
pi = P1(i, :);
pj = P1(j, :);
while (diff > epsilon)
    
    scra1 = mu2 * W';
    scra2 = mu1 * W;
    
    new1 = pi./scra1;
    new1 = new1/sum(new1);
    
    new2 = pj./scra2;
    new2 = new2/sum(new2);
    
    diff = max(max(abs(new1-mu1), abs(new2-mu2)));
    
    mu1 = new1;
    mu2 = new2;
    
end
end

function DI = compute_di(i, j, W, mu1, mu2, Pia)
tiny = 1.0e-100;
Pdir = W.*(mu1'*mu2);
Pdir = Pdir / sum(sum(Pdir));
Pfac = Pia(i, :)' * Pia(j, :);
DI = trace(Pdir' * log((Pdir+tiny)./(Pfac+tiny)));
end

function [index_of_interest, range_start, range_end] = find_seq_of_interest(full_alignment, seqid_of_interest)
index_of_interest = -1;
for scan_index = 1:size(full_alignment,1)
    [seqid, range_start, range_end] = split_uniprot_id(full_alignment(scan_index).Header);
    if (strcmp(seqid,seqid_of_interest) == 1)
        index_of_interest = scan_index;
        break
    end
end
if (index_of_interest == -1)
    error(['Error: could not find sequence of interest (' seqid_of_interest ') in multiple sequence alignment']);
end
end

function [seqid, range_start, range_end] = split_uniprot_id(pfam_uniprot_range_line)
slashposition = findstr('/', pfam_uniprot_range_line);
if (size(slashposition,2) ~= 1 || slashposition == 1 || slashposition == size(pfam_uniprot_range_line,2))
    error(['Error: could not parse (slash error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
seqid = pfam_uniprot_range_line(1:slashposition - 1);
range_start=1;
range_end=[];
%rangestring = pfam_uniprot_range_line(slashposition + 1:size(pfam_uniprot_range_line,2));
%hyphenposition = findstr('-', rangestring);
%if (size(hyphenposition,2) ~= 1 || hyphenposition == 1 || hyphenposition == size(rangestring,2))
%    error(['Error: could not parse (hyphen error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
%end
%range_start = str2num(rangestring(1:hyphenposition - 1));
%range_end = str2num(rangestring(hyphenposition + 1 : size(rangestring,2)));
%if (isempty(range_start) || isempty(range_end))
%    error(['Error: could not parse (range start/end) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
%end
end

function [focuscolumnlist, uniprotoffsetlist] = scan_sequence_of_interest_for_focus_columns(sequence_of_interest, range_of_interest_start, letter2number_map)
focuscolumnlist = [];
uniprotoffsetlist = [];
next_uniprotoffset = range_of_interest_start;
for pos=1:size(sequence_of_interest,2)
    residuecode = letter2number_map(sequence_of_interest(pos));
    if (residuecode == 0)
        error(['Error: sequence of interest contains undefined residues:' sequence_of_interest]);
    end
    if (residuecode == -1)
        error(['Error: sequence of interest contains ambiguous residues:' sequence_of_interest]);
    end
    if (residuecode > 1)
        focuscolumnlist = [focuscolumnlist pos];
        uniprotoffsetlist = [uniprotoffsetlist next_uniprotoffset];
    end
    if (residuecode == -2 || residuecode > 1)
        next_uniprotoffset = next_uniprotoffset + 1;
    end
end
end

function letter2number_map = create_letter2number_map()
letter2number_map(256) = 0; %initiallize all bytes to 0
letter2number_map('-') = 1;
letter2number_map('A') = 2;
letter2number_map('C') = 3;
letter2number_map('D') = 4;
letter2number_map('E') = 5;
letter2number_map('F') = 6;
letter2number_map('G') = 7;
letter2number_map('H') = 8;
letter2number_map('I') = 9;
letter2number_map('K') = 10;
letter2number_map('L') = 11;
letter2number_map('M') = 12;
letter2number_map('N') = 13;
letter2number_map('P') = 14;
letter2number_map('Q') = 15;
letter2number_map('R') = 16;
letter2number_map('S') = 17;
letter2number_map('T') = 18;
letter2number_map('V') = 19;
letter2number_map('W') = 20;
letter2number_map('Y') = 21;
letter2number_map('B') = -1; %ambiguous : skip sequences containing these
letter2number_map('Z') = -1; %ambiguous : skip sequences containing these
letter2number_map('J') = -1; %ambiguous : skip sequences containing these
letter2number_map('X') = -1; %ambiguous : skip sequences containing these
letter2number_map('U') = -1; %non-standard : skip sequences containing these
letter2number_map('O') = -1; %non-standard : skip sequences containing these
letter2number_map('a') = -2; %non-conserved: skip in seq of interest
letter2number_map('c') = -2; %non-conserved: skip in seq of interest
letter2number_map('d') = -2; %non-conserved: skip in seq of interest
letter2number_map('e') = -2; %non-conserved: skip in seq of interest
letter2number_map('f') = -2; %non-conserved: skip in seq of interest
letter2number_map('g') = -2; %non-conserved: skip in seq of interest
letter2number_map('h') = -2; %non-conserved: skip in seq of interest
letter2number_map('i') = -2; %non-conserved: skip in seq of interest
letter2number_map('k') = -2; %non-conserved: skip in seq of interest
letter2number_map('l') = -2; %non-conserved: skip in seq of interest
letter2number_map('m') = -2; %non-conserved: skip in seq of interest
letter2number_map('n') = -2; %non-conserved: skip in seq of interest
letter2number_map('p') = -2; %non-conserved: skip in seq of interest
letter2number_map('q') = -2; %non-conserved: skip in seq of interest
letter2number_map('r') = -2; %non-conserved: skip in seq of interest
letter2number_map('s') = -2; %non-conserved: skip in seq of interest
letter2number_map('t') = -2; %non-conserved: skip in seq of interest
letter2number_map('v') = -2; %non-conserved: skip in seq of interest
letter2number_map('w') = -2; %non-conserved: skip in seq of interest
letter2number_map('y') = -2; %non-conserved: skip in seq of interest
letter2number_map('b') = -2; %non-conserved: skip in seq of interest
letter2number_map('z') = -2; %non-conserved: skip in seq of interest
letter2number_map('j') = -2; %non-conserved: skip in seq of interest
letter2number_map('x') = -2; %non-conserved: skip in seq of interest
letter2number_map('u') = -2; %non-conserved: skip in seq of interest
letter2number_map('o') = -2; %non-conserved: skip in seq of interest
letter2number_map('.') = -3; %non-conserved: skip in seq of interest, do not advance position
end

function number2letter_map = create_number2letter_map()
number2letter_map(1) = '-';
number2letter_map(2) = 'A';
number2letter_map(3) = 'C';
number2letter_map(4) = 'D';
number2letter_map(5) = 'E';
number2letter_map(6) = 'F';
number2letter_map(7) = 'G';
number2letter_map(8) = 'H';
number2letter_map(9) = 'I';
number2letter_map(10) = 'K';
number2letter_map(11) = 'L';
number2letter_map(12) = 'M';
number2letter_map(13) = 'N';
number2letter_map(14) = 'P';
number2letter_map(15) = 'Q';
number2letter_map(16) = 'R';
number2letter_map(17) = 'S';
number2letter_map(18) = 'T';
number2letter_map(19) = 'V';
number2letter_map(20) = 'W';
number2letter_map(21) = 'Y';
end

function [data, seq] = myfastaread(filename,varargin)
%FASTAREAD reads FASTA format file.
%
%   S = FASTAREAD(FILENAME) reads a FASTA format file FILENAME, returning
%   the data in the file as a structure. FILENAME can also be a URL or
%   MATLAB character array that contains the text of a FASTA format file.
%   S.Header is the header information. S.Sequence is the sequence stored
%   as a string of characters.
%
%   [HEADER, SEQ] = FASTAREAD(FILENAME) reads the file into separate
%   variables HEADER and SEQ. If the file contains more than one sequence,
%   then HEADER and SEQ are cell arrays of header and sequence information.
%
%   FASTAREAD(...,'IGNOREGAPS',TF) removes any gap symbol ('-' or '.')
%   from the sequence(s) when TF is true. Default is false.
%
%   FASTAREAD(...,'BLOCKREAD', M) allows you to read in a single entry or
%   block of entries from a file containing multiple sequences. If M is a
%   scalar then the M'th entry in the file is read. If M is a two element
%   vector then the block of entries starting at entry M(1) and ending at
%   entry M(2) will be read.  Use Inf for M(2) to read all entries in the
%   file starting at position M(1).
%
%   FASTAREAD(...,'TRIMHEADERS',TF) trims the header after the first
%   whitespace when TF is true. White space characters include a horizontal
%   tab (char(9)) and a space (char(32)). Default is false.
%
%   FASTA format specified here:
%   http://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml
%
%   Examples:
%
%       % Read the sequence for the human p53 tumor gene.
%       p53nt = fastaread('p53nt.txt')
%
%       % Read the sequence for the human p53 tumor protein.
%       p53aa = fastaread('p53aa.txt')
%
%       % Read a block of entries from a file
%       pf2_5_10 = fastaread('pf00002.fa','blockread',[ 5 10], ...
%                            'ignoregaps',true)
%
%   See also EMBLREAD, FASTAINFO, FASTAWRITE, FASTQINFO, FASTQREAD,
%   FASTQWRITE, GENBANKREAD, GENPEPTREAD, HMMPROFDEMO, MULTIALIGNREAD,
%   MULTIALIGNVIEWER, SEQPROFILE, SEQTOOL, SFFINFO, SFFREAD.

%   Copyright 2002-2011 The MathWorks, Inc.
%   $Revision: 1.15.4.23 $  $Date: 2011/02/09 18:49:21 $

% check input is char
% in a future version we may accept also cells
if ~ischar(filename)
    error(message('bioinfo:fastaread:InvalidInput'))
end

% default
ignoreGaps = false;
trimh = false;
blockRead = false;
% get input arguments
if  nargin > 1
    if rem(nargin,2) ~= 1
        error(message('bioinfo:fastaread:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'ignoregaps','blockread','trimheaders'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:fastaread:UnknownParameterName', pname));
            %elseif length(k)>1
            %    error('bioinfo:fastaread:AmbiguousParameterName',...
            %        'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % ignore gaps
                    ignoreGaps = opttf(pval);
                    if isempty(ignoreGaps)
                        error(message('bioinfo:fastaread:IgnoreGapsNotLogical', upper( char( okargs( k ) ) )));
                    end
                case 2  % range
                    range = pval;
                    if ~isnumeric(range) || numel(range)> 2 || isempty(range)
                        error(message('bioinfo:fastaread:BadBlockRange'))
                    end
                    blockRead = true;
                    range = sort(range);
                case 3 % trimheaders
                    trimh = opttf(pval, okargs{k}, mfilename);
            end
        end
    end
end


if size(filename,1)>1  % is padded string
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    ftext = cell(size(filename,1),1);
    for i=1:size(filename,1)
        ftext(i,1)=strtrim(strread(filename(i,:),'%s','whitespace','','delimiter','\n'));
    end
    % try then if it is an url
elseif (strfind(filename(1:min(10,end)), '://'))
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    if (~usejava('jvm'))
        error(message('bioinfo:fastaread:NoJava'))
    end
    try
        ftext = urlread(filename);
    catch allExceptions
        error(message('bioinfo:fastaread:CannotReadURL', filename));
    end
    ftext = strread(ftext,'%s','delimiter','\n');
    
    % try then if it is a valid filename
elseif  (exist(filename,'file') || exist(fullfile(pwd,filename),'file'))
    if blockRead
        blockText = getfileblock(filename,range,'>');
        try
            ftext = strread(blockText,'%s','delimiter','\n');
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:fastaread:BlockTooBig'));
            else
                rethrow(theErr);
            end
        end
    else %=== read entire file
        fid = fopen(filename);
        c = onCleanup(@()fclose(fid));
        try
            ftext = textscan(fid,'%s','delimiter','\n');
            ftext = ftext{:};
        catch theErr
            if strcmpi(theErr.identifier,'MATLAB:nomem')
                error(message('bioinfo:fastaread:FileTooBig'));
            else
                rethrow(theErr);
            end
        end
    end
else  % must be a string with '\n', convert to cell
    if blockRead
        warning(message('bioinfo:fastaread:IgnoredRange'))
    end
    ftext = strread(filename,'%s','delimiter','\n');
end

% it is possible that there will be multiple sequences
commentLines = strncmp(ftext,'>',1);

if ~any(commentLines)
    error(message('bioinfo:fastaread:FastaNotValid'))
end

numSeqs = sum(commentLines);
seqStarts = [find(commentLines); size(ftext,1)+1];
data(numSeqs,1).Header = '';

try
    for theSeq = 1:numSeqs
        % Check for > symbol ?
        data(theSeq).Header = ftext{seqStarts(theSeq)}(2:end);
        firstRow = seqStarts(theSeq)+1;
        lastRow = seqStarts(theSeq+1)-1;
        numChars = cellfun('length',ftext(firstRow:lastRow));
        numSymbols = sum(numChars);
        data(theSeq).Sequence = repmat(' ',1,numSymbols);
        pos = 1;
        for i=firstRow:lastRow,
            str = strtrim(ftext{i});
            len =  length(str);
            if len == 0
                break
            end
            data(theSeq).Sequence(pos:pos+len-1) = str;
            pos = pos+len;
        end
        data(theSeq).Sequence = strtrim(data(theSeq).Sequence);
        if ignoreGaps
            data(theSeq).Sequence = strrep(data(theSeq).Sequence,'-','');
            data(theSeq).Sequence = strrep(data(theSeq).Sequence,'.','');
        end
    end
    
    % trim headers
    if trimh
        for i = 1:numSeqs
            data(i).Header = sscanf(data(i).Header,'%s',1);
        end
    end
    
    % in case of two outputs
    if nargout == 2
        if numSeqs == 1
            seq = data.Sequence;
            data = data.Header;
        else
            seq = {data(:).Sequence};
            data = {data(:).Header};
        end
    end
    
catch allExceptions
    error(message('bioinfo:fastaread:IncorrectDataFormat'))
end
end
