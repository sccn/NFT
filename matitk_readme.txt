Setup for Windows:
==================
1. copy the DLLs to a desired location
2. set search path of matlab or change current directory to the location of the dll

Example MATLAB commands:
========================
matitk('?')
matitk('f')
matitk('s')
matitk('r')
%------
load mri;
D=squeeze(D);
b=matitk('FCA',[5 0.0625 3],double(D));
c=matitk('SCC',[1.4 10 255],double(b),[],[102 82 25]);subplot(131);
imagesc(squeeze(D(:,:,15)));
axis image;
colormap gray
subplot(132);
imagesc(squeeze(b(:,:,15)));
axis image;
colormap gray
subplot(133);
imagesc(squeeze(c(:,:,15)));
axis image;
colormap gray
%------

A detail usage guide can be found at:
http://www.sfu.ca/~vwchu/matitkusage.html

More Information
================

Please report any comments/bugs to: Vincent Chu <vwchu@sfu.ca>, Ghassan Hamarneh <hamarneh@cs.sfu.ca>

Additional details about the project can be found at:
http://mial.fas.sfu.ca/researchProject.php?s=308


