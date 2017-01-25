%
%  main.m
%  SfSM: High Quality Structure from Small Motion for Rolling Shutter Cameras
% 
%
%  Created by Sunghoon Im on 2017. 1. 25..
%  Copyright @ 2017 Sunghoon Im. All rights reserved.
%

clc; clear; close all;

addpath('harris');
setting.srcdir='60D';
setting.type_data='MOV';
setting.alpha=0.5;
setting.InitDepth = 100;

flist=dir(fullfile(setting.srcdir,['*.' setting.type_data]));
for i=1:length(flist)
    setting.cname=flist(i).name(1:end-4);
    A = SfSM(setting);
    A.SMBARS();
end

            

