function [rec, det, lam, corm, Rec_matrix] = Trqa(xcoor, ycoor, dur, Rshow, radius, linelength)
% TRQA This script calculates four recurrence analysis measures based on fixation locations and duration.
%   To run this code and obtain recurrence measures data should be provided.
%   If using the code, please cite the paper 'Recurrence quantification analysis reveals eye-movement behavior differences between experts and novices
%   P Vaidyanathan, J Pelz, C Alm, P Shi, A Haake - Proceedings of the Symposium on Eye Tracking Research Application, ETRA 2014

%   This function calculates the recurrence, determinism, laminarity and
%   center of recurrence mass as defined by the equations in Anderson et.
%   al.,(Behavioral Res, 2013).
%  INPUT:   xcoor      --- x coordinates for eye movement fixations (column vector Nx1)
%           ycoor      --- y coordinates for eye movement fixations (column vector Nx1)
%           dur        --- fixation duration (column vector - Nx1)
%           Rshow      --- whether to plot the recurrence values 
%           radius     --- dist. b/w fixations for recurrence (default = 64, based on foveal approximation)
%           linelength --- line length for vert., horiz. and diagnoal lines (scalar, default = 2)
%  OUTPUT:  rec        --- recurrence measure (scalar)
%           det        --- determinism measure (scalar)
%           lam        --- laminarity measure (scalar)
%           corm       --- center of recurrence mass (scalar)
%           REC_matrix --- recurrence matrix for plotting (matrix)

% Default parameters
if nargin <4
     Rshow = 0;
end
if nargin <5
    radius = 64;       
end
if nargin < 6
    linelength = 2;   
end



% Number of fixations
NofFix = size(xcoor,1);

% Distance matrix - distance between each fixation to the other fixation
dist_matrix = sqrt((bsxfun(@minus,xcoor,xcoor')).^2+(bsxfun(@minus,ycoor,ycoor')).^2);

% sum of fixation durations for two fixations
SumofDur  = meshgrid(dur)+meshgrid(dur)';
% Threshold for the given radius (fixation distance)
Rec_matrix = dist_matrix <= radius;
Rec_matrix_T = Rec_matrix.*SumofDur;

% Get only the upper triangle without the lin of incidence
UT_Rec_matrix = triu(Rec_matrix,1);%binary
UT_Rec_matrix_T = triu(Rec_matrix_T,1);%with fixation durations

% Recurrence measures
R_T = sum(UT_Rec_matrix_T(:));

% total fixation duration
T = sum(dur);

% Diagonal lines and sum the fixation duration between fixations that are
% in a diagonal
[R,I] = spdiags(UT_Rec_matrix);
[Rdiagonals,Idiagonals] = spdiags(UT_Rec_matrix_T);
Diag_Lines = [];
Diag_Values = [];
for i = 1:length(I)
    d = [0;R(:,i);0];
    diagchange = diff(d);
    dl = find(diagchange==-1) - find(diagchange==1);
    [p,~] = find(diagchange==-1); [q,~] = find(diagchange==1);
    p = p.*(dl>=linelength); q = q.*(dl>=linelength);
    c = [nonzeros(q) nonzeros(p)-1];
    for ii = 1:size(c,1)
       Diag_Values =  [Diag_Values;Rdiagonals(c(ii,1):c(ii,2),i)];
    end
    dl = dl(dl >= linelength);
    Diag_Lines = [Diag_Lines;dl];
end
Diag_Lines = sum(Diag_Lines);
Diag_Values = sum(Diag_Values);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% Horizontal & Vertical lines
VertHorz_Lines = [];
VertHorz_Values = [];
% remove line of incidence - main diagonal but keep rest
paddedRM = [zeros(NofFix,1) Rec_matrix-eye(NofFix) zeros(NofFix,1)];
paddedRM_T = Rec_matrix_T-eye(NofFix);
paddedRM=max(paddedRM,0);
for i = 1:NofFix
    v = paddedRM(i,:);
    vchange = diff(v);
    vl = find(vchange==-1)-find(vchange==1);
    [~,p] = find(vchange==-1); [~,q] = find(vchange==1);
    p = p.*(vl>=linelength); q = q.*(vl>=linelength);
    c = [nonzeros(q) nonzeros(p)-1];
    for ii = 1:size(c,1)
       VertHorz_Values =  [VertHorz_Values paddedRM_T(i,c(ii,1):c(ii,2))];
    end
    vl = vl(vl >= linelength);
    VertHorz_Lines = [VertHorz_Lines vl]; 
end
VertHorz_Lines = sum(VertHorz_Lines);
VertHorz_Values = sum(VertHorz_Values);

% Recurrence REC
rec = 100*R_T/((NofFix-1)*T);

% Determinism DET
det = 100*Diag_Values/R_T;

% Laminarity LAM
lam = 100*VertHorz_Values/(2*R_T);

% CORM - center of recurrence mass
corm = NaN;
if ~isempty(Rdiagonals)&& ~isempty(Idiagonals)
corm = 100*sum(sum(Rdiagonals)'.*Idiagonals)/((NofFix-1)^2*R_T);
end

% show the recurrence plot
if Rshow == 1
    image(Rec_matrix_T*255);
    colormap(gray),set(gca,'YDir','normal');
end    
end



