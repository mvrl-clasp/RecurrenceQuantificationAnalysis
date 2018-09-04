% This is a demo of the scripts rqa.m and Trqa.m to perform recurrence quantification
% analysis comparison between experts and novices. The output values would
% be: rec = 3.5/3.4, det = 14.6/12.12, lam = 33.7/30, corm = 29.1/0.3 and the output rqa
% plot would look like DemoFigure.fig

% Read DemoFixation.csv
% FixationData =
% csvread('~/desktop/vlsa/vlsadata/clusters/MSFC/fixationn.csv');
FixationData = csvread('DemoFixation.csv');
xcoor = FixationData(:,1);
ycoor = FixationData(:,2);
dur   = FixationData(:,3);

% Define your parameters
Rshow      = 1;
radius     = 64;
linelength = 2;

% Uncomment the following line to run rqa
%[rec, det, lam, corm, Rec_matrix] = rqa(xcoor, ycoor, Rshow, radius, linelength);

% Uncomment the following line to run rqa with fixation duration
 [rec, det, lam, corm, Rec_matrix] = Trqa(xcoor, ycoor, dur, Rshow, radius, linelength);