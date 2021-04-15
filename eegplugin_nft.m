% eegplugin_nft() - NFT plugin 
%
% Usage:
%   >> eegplugin_nft(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Author: Zeynep Akalin Acar, SCCN, INC, UCSD

% Copyright (C) 2009-, Zeynep Akalin Acar
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

% $Log: eegplugin_nft.m,v $

function vers = eegplugin_nft(fig, trystrs, catchstrs)
    
    
    vers = 'nft2.3';
    if nargin < 3
        error('eegplugin_nft requires 3 arguments');
    end;
    
    % add NFT folder to path
    % -----------------------
    if ~exist('eegplugin_nft')
        p = which('eegplugin_nft');
        p = p(1:findstr(p,'eegplugin_nft.m')-1);
        addpath(p);
    end;

    if ~exist('NFT')
    	disp('Warning: Please install NFT and make sure it is added to MATLAB path');
    end;

    % find tools menu
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 
    % tag can be 
    % 'import data'  -> File > import data menu
    % 'import epoch' -> File > import epoch menu
    % 'import event' -> File > import event menu
    % 'export'       -> File > export
    % 'tools'        -> tools menu
    % 'plot'         -> plot menu
    
    % menu callback commands
    % ----------------------
    comrun1   = [ 'NFT(''EEGstruct'',EEG)'];
    %comrun2   = [ 'xxxxxxxxxxxxxx put your command in here xxxxxxxxxxxxxxx'  ];
    
    % create menus
    % ------------
    submenu = uimenu( menu, 'Label', 'NFT plugin', 'separator', 'on');
    uimenu( submenu, 'Label', 'Start NFT GUI', 'CallBack', comrun1);
    %uimenu( submenu, 'Label', 'Inverse Problem with NFT', 'CallBack', comrun2);
