#!/usr/bin/octave

addpath('../utils')
close all; clc;

%Nheader=1; % Number of header lines in timeprogress files
%ycol=4;    %
if  nargin==0
  printf("Usage: ./get_StickyDiffusivity.m <Ze> <Zs> <p> <tauS>\n");
else

arg_list = argv();
ZE  =str2num(arg_list{1}); % Number of entanglements per chain
ZS  =str2num(arg_list{2}); % Number of stickers per chain
p   =str2num(arg_list{3}); % Fraction of bound stickers
tauS=str2num(arg_list{4}); % Sticker lifetime

DR=1.0/(3*ZE*pi^2);
DSR=StickyDiffusionCoefficient(p, ZE, ZS, tauS); # Units of DR
fprintf( '%e\n', DSR*DR);
end

