% main script
function [B, B0, Ic, ImaskClean] = segm_inhomogeneity_correction(I0)
%I0 = imread([filepath filename]);
%I0=imread('sl_97.tiff');
I0 = single(I0);
Imin = min(I0(:));
Imax = max(I0(:));
I0 = I0 - Imin;
Ii = I0 / (Imax-Imin)*255;


% 3. Thresholding
% REGION GROWING FROM CORNER
gain = 3;
[bck_mask,bck_mean,bck_std] = NeckBackground(Ii,0,gain);
deg = sum(sum(bck_mask))/size(bck_mask,1)/size(bck_mask,2);
while deg<0.2 & gain<6
    gain = gain+1;
    [bck_mask,bck_mean,bck_std] = NeckBackground(Ii,0,gain);
    deg = sum(sum(bck_mask))/size(bck_mask,1)/size(bck_mask,2);
end
pa = 2;
pb = 4;
fuzzy = [pa  pb];
% disp(['     Fuzzy rule: ' num2str(fuzzy)])
I = medfilt2(Ii, [5 5]);
Imask = segm_inhomog_neckmask(I, bck_mean + 0.655 * bck_std * fuzzy);
ImaskClean = NeckMaskCLean( Imask ); % !!!!!

% 4. Initialize bias field    
B0_order = 3; 
B0 = PolyMaskFilter(double(Ii), B0_order, ImaskClean);
sig = 31;
G = fspecial('gaussian', 3*sig+1, sig);
temp = imfilter(double(Ii), G, 'same', 'replicate');
data_B0 = max(0.5 * temp, B0);

% 5. Filter
iterations = 9;
sigma = 10; %?
If = anisoOS(Ii, 'tukeyPsi', sqrt(2) * sigma, iterations, 1, data_B0);

% 6. Correct Intensity 
options = BiasCorrLEMSS2D;
a = 20;
b = 20;
options.Nknots = [a b];
options.NiterMax = 1; %4;
options.Bgain = 0.8;
options.flag_display = 0;

algo = 'LEMS2D';
  
B = BiasCorrLEMSS2D(If, Imask.*ImaskClean, bck_mean, bck_std, options, data_B0);
Ic = Ii./B.*Imask + (1-Imask).*Ii;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Imask4 = NeckMaskCLean(Imask);
% NeckMaskCLean:    clean a neck mask by removing small elements and opening
%
%   Imask4 = NeckMaskCLean(Imask);
%
% OS CWRU, 19-jan-04


% --- remove small stuff.
Imask = imopen(Imask,strel('disk',1,0));

% --- keep the biggest CCA
Imask2 = bwlabel(Imask,8);
s  = regionprops(Imask2, 'Area');
area = [s.Area];
[dummy,idx] = max(area);
Imask3 = (Imask2==idx);

% --- close the holes
Imask4 = imclose(Imask3,strel('disk',2,0));
Imask4 = bwmorph(Imask4,'fill');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask,V_mean,V_std,Vfad] = NeckBackground(V,flagdisplay,gain);
% NeckBackground:       background from neck image by region growing from
%   [mask,V_mean,V_std,Vfad] = NeckBackground(V,flagdisplay,gain);
%   gain: multiplicative gain of STD for threshold (default=3)
% 
% OS CWRU, 05-jun-03

if ~exist('flagdisplay'),
    flagdisplay=0;
end
if ~exist('gain'),
    gain = 3;
end



if flagdisplay,
%     colormap(gray)
%     subplot(121)
%     imagesc(V),axis image, axis off
%     subplot(122)
    disp('--- Starting region growing from corner')
end

% clear all
SE3 = strel('disk',3,0);
SE1 = strel('disk',1,0);

threshold = 0.001;


[nr,nc,nf] = size(V);

% --- filter image using anisotropic diffusion
Vfad = anisoOS(V,'tukeyPsi',0.5*mean2(V),30);    

% --- get a seed area
Vgrow = im2bw(0*V);
Vgrow(1,1) = 1;
Vgrow(1,end) = 1;
new_voxel = 1000;

% --- main loop until no more voxel added
goahead = 1;
flipflop = 1;
idx = 1;
% if flagdisplay,
%     disp('  Start the region growing')
% end
Vinit = 1e6*ones(size(V));
while goahead,
    idx = idx+1;
    
    % --- compute the statistics of the area
    %     V_mean = sum(sum(Vfad.*Vgrow)) / sum(sum(Vgrow));
    %     V_std = sqrt(sum(sum((Vfad-V_mean).^2.*Vgrow)))/...
    %         (sum(sum(Vgrow))-1);
    V_mean = mean(Vfad(Vgrow));
    V_std = std(Vfad(Vgrow));
    
    threshold = V_mean + gain* V_std;
    %     threshold = 2* V_mean;
    
    % --- ... and from within a slice
    if new_voxel>30,
        SE = SE3;
    else
        SE = SE1;
    end
    Vnewslice = imdilate(Vgrow,SE);
    
    % --- The new label voxel to be tested
    Vnew = im2bw(Vnewslice.*(1-Vgrow) );
    
    % --- criteria for each new voxel
    Vcrit = Vinit;
    Vtemp = (Vnew.*Vfad - V_mean).^2;
    Vcrit(Vnew) = Vtemp(Vnew);
    
    % --- get the new good voxels
    Vgood = 0*V;
    Vgood = Vcrit < threshold;
    
    % --- grow the region
    Vgrow = or(Vgrow , Vgood);
    
    % --- total number of voxel in the regions
    NbVox = sum(sum(Vgrow));
    
    % --- fill the holes
    Vgrow = bwmorph(Vgrow,'fill');  % remove holes
    
    % --- continue of the region has grown
    new_voxel = sum(sum(Vgood));
    goahead = new_voxel > 0;
    
    if flagdisplay>1,
        disp([' Mean:' num2str(V_mean) '  Number of new vox:' num2str(new_voxel)...
                '  STD:' num2str(V_std) '  Threshold:' num2str(threshold) ])
        imagesc(Vgrow),axis image, axis off
        drawnow
    end
end
% --- close the small holes
mask = imclose(Vgrow,strel('disk',4,0));
V_mean = mean(V(mask));
V_std = std(V(mask));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    [Imask,Ibck] = segm_inhomog_neckmask(I,alpha,corner,percent);
% NeckMask:     get a mask for the background 
%
%   [Imask,Ibck] = NeckMask(I,alpha,corner,percent);
%       alpha: threshold for masking 0.1 <=> 10 % upper histogramme is one
%               0, uses 10% of the corner
%           [amin amax] : compute weighted mask with
%                       Imask(i,j) = 0 if I<amin
%                       Imask(i,j) = 1 if I>amax   
%                       linear interpolation in between
%       corner: 1 uses top corners with percent of the dimension
%               2 uses 4 corners with percent of the dimension
%               


I   = single(I);
flagW = 0;
% --- default parameters
if ~exist('alpha'),
    alpha = 0;
end
if length(alpha) > 1,
    flagW = 1;
    corner = 0;
    percent = 0;
end
if ~exist('corner'),
    corner = 2;
    percent = 0.2;
end
if ~exist('percent'),
    percent = 0.2;
end

if ~flagW,
    if alpha==0,
        % --- use the corner to get the background
        [nr,nc] = size(I);
        r_cner  = round(nr*percent);    % 10% of the image size
        c_cner  = round(nc*percent);
        Ic1     = I(1:r_cner,           1:c_cner);
        Ic2     = I(1:r_cner,           end-c_cner:end);
        Ic3     = I(end-r_cner:end,     1:c_cner);
        Ic4     = I(end-r_cner:end,     end-c_cner:end);
        
        Ibck    = [Ic1];
        if corner==4,
            %         disp('four corners')
            Ibck    = [Ic1 Ic2;Ic3 Ic4];
        elseif corner==2,
            %         disp('2 top corners')
            Ibck    = [Ic1 Ic2];
        end
        
        Im      = mean2(Ibck);
        Is      = std2(Ibck);
        [wx,wy] = find( I > (Im+3*Is));
        idx     = sub2ind(size(I),wx,wy);
        
    else,
        
        Imax = max(max(I));
        Imin = min(min(I));
        threshold   = Imin + (Imax-Imin)*alpha;
        
        [wx,wy] = find( I > threshold);
        idx     = sub2ind(size(I),wx,wy);
        Ibck    = [];
    end
else,
    % --- Weighted mask
    if alpha(1)<1,  % assume it is a fraction of the histogram
        Imax = max(max(I));
        Imin = min(min(I));
        Idelta = Imax-Imin;
        Itl = Imin + alpha(1)*Idelta;
        Ith = Imin + alpha(2)*Idelta;
    else, %assume the thresholds are in absolute value
        Itl = alpha(1);
        Ith = alpha(2);
    end
    
    
    Imask = 0*I;    % below amin
    Imask(I>=Ith) = 1;  % above amax
    
    inbetween = 0*I;
    inbetween = (I>Itl) & (I<Ith);
    Imask(inbetween) = ( I(inbetween)-Itl )/ (Ith-Itl);
end

Ibck = I.*(1-Imask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = anisoOS(input,psiFunction,sigma,iterations,lambda,B)
%
% result = aniso(input,[psiFunction],[sigma],[iterations],[lambda])
%
% Anisotropic diffusion, following the robust statistics
% formulation laid out in Black et al, IEEE Trans on Im Proc,
% 7:421-432, 1998.
%
% Examples of how to run this function are included at the end of
% aniso.m 
%
%   input: input image
%   psiFunction: influence function that determines how the
%      diffusion depends on the local image gradient.  Three
%      example psi functions are:
%         linearPsi.m
%         lorentzianPsi.m.
%         tukeyPsi.m
%      Default is tukeyPsi, that makes this the same as one of
%      original two algorithms proposed by Perona et al.  But
%      tukeyPsi is often a better choice because it does a better
%      job of maintaining the sharpness of an edge.  LinearPsi
%      gives standard linear diffusion, i.e., shift-invariant
%      convolution by a Gaussian. 
%   sigma: scale parameter on the psiFunction.  Choose this
%      number to be bigger than the noise but small than the real
%      discontinuties. [Default = 1]
%   iterations: number of iterations [Default = 10]
%   lambda: rate parameter [Default = 1/4] To approximage a
%      continuous-time PDE, make lambda small and increase the
%      number of iterations.
%
% DJH, 2/2000
%
% modified by OS 28sep03 to removed gradient from B.


if ~exist('psiFunction','var')
   psiFunction = 'tukeyPsi';
end
if ~exist('sigma','var')
   sigma = 1;
end
if ~exist('iterations','var')
   iterations = 10;
end
if ~exist('lambda','var')
   lambda = 0.25;
end
lambda = lambda/4;

% Initialize result
result = input;

% Indices for the center pixel and the 4 nearest neighbors
% (north, south, east, west)
[m,n] = size(input);
rowC = [1:m];         rowN = [1 1:m-1];    rowS = [2:m m];
colC = [1:n];         colE = [1 1:n-1];    colW = [2:n n];

% --- take into account B if it exists
if exist('B','var')
   northB = B(rowN,colC)-B(rowC,colC);
   southB = B(rowS,colC)-B(rowC,colC);
   eastB = B(rowC,colE)-B(rowC,colC);
   westB = B(rowC,colW)-B(rowC,colC);
end



for i = 1:iterations
   % Compute difference between center pixel and each of the 4
   % nearest neighbors.
   north = result(rowN,colC)-result(rowC,colC);
   south = result(rowS,colC)-result(rowC,colC);
   east = result(rowC,colE)-result(rowC,colC);
   west = result(rowC,colW)-result(rowC,colC);
   if exist('B','var')
       north = north - northB;
       south = south - southB;
       east = east - eastB;
       west = west - westB;
   end
       
   % Evaluate the psiFunction for each of the neighbor
   % differences and add them together.  If the local gradient is
   % small, then the psiFunction should increase roughly linearly
   % with the neighbor difference.  If the local gradient is large
   % then the psiFunction should be zero (or close to zero) so
   % that large gradients are ignored/treated as outliers/stop the
   % diffusion.
   psi = eval([psiFunction,'(north,sigma)']) + ...
      eval([psiFunction,'(south,sigma)']) + ...
      eval([psiFunction,'(east,sigma)']) + ...
      eval([psiFunction,'(west,sigma)']);
   % Update result
   result = result + lambda * psi;
end;
return


function y = linearPsi(x,sigma)
y = 2*x;
return;


function y = lorentzianPsi(x,sigma)
y = (2*x)./(2*sigma^2 + abs(x).^2);
return


function y = tukeyPsi(x,sigma)
y = zeros(size(x));
id = (x > -sigma) & (x < sigma);
xid = x(id);
y(id) = xid.*((1-(xid/sigma).^2).^2);
return


% Test, debug, examples

% Noisy step
step = 1+[ones([32 64]) ; zeros([32 64])];
noise = 0.2 * randn(size(step));
im = (step + noise);

% Small number of iterations.  Not much difference.
resultLin10 = aniso(im,'linearPsi',0,10);
resultLor10 = aniso(im,'lorentzianPsi',0.5,10);
resultTuk10 = aniso(im,'tukeyPsi',0.5,10);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor10,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk10,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin10(:,32) resultLor10(:,32) resultTuk10(:,32)]);

% More iterations.  Note that tukeyPsi is much more robust
% because it is a fully redescending influence function (see
% Black et al).  The upshot of this is that it maintains the
% sharpness of the edge even after a large number of iterations
% whereas lorentzianPsi does not.
resultLin100 = aniso(im,'linearPsi',0,100);
resultLor100 = aniso(im,'lorentzianPsi',0.5,100);
resultTuk100 = aniso(im,'tukeyPsi',0.5,100);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor100,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk100,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin100(:,32) resultLor100(:,32) resultTuk100(:,32)]);

% Lots o' iterations.  Tukey converges wheres lorentzian and
% linear continue to blur more and more.
resultLin1000 = aniso(im,'linearPsi',0,1000);
resultLor1000 = aniso(im,'lorentzianPsi',0.5,1000);
resultTuk1000 = aniso(im,'tukeyPsi',0.5,1000);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor1000,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk1000,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin1000(:,32) resultLor1000(:,32) resultTuk1000(:,32)]);

% Tradeoff between lambda and iterations.  All of these should
% give the same result
result5 = aniso(im,'tukeyPsi',0.5,5,0.5);
result10 = aniso(im,'tukeyPsi',0.5,10,0.25);
result100 = aniso(im,'tukeyPsi',0.5,100,0.025);
result1000 = aniso(im,'tukeyPsi',0.5,1000,0.0025);
figure(1); clf;
plot([im(:,32) result5(:,32) result10(:,32) result100(:,32) result1000(:,32)]);


% ------------- test bias field effect -----------------------
% Noisy step
step = 1+[ones([32 64]) ; zeros([32 64])];
[x,y] = meshgrid(1:64,1:64);
B = cos(x/3);
noise = 0.2 * randn(size(step));
im = (step + noise);

% Small number of iterations.  Not much difference.
resultLin10 = aniso(im,'linearPsi',0,10);
resultLor10 = aniso(im,'lorentzianPsi',0.5,10);
resultTuk10 = aniso(im,'tukeyPsi',0.5,10);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor10,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk10,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin10(:,32) resultLor10(:,32) resultTuk10(:,32)]);

% More iterations.  Note that tukeyPsi is much more robust
% because it is a fully redescending influence function (see
% Black et al).  The upshot of this is that it maintains the
% sharpness of the edge even after a large number of iterations
% whereas lorentzianPsi does not.
resultLin100 = aniso(im,'linearPsi',0,100);
resultLor100 = aniso(im,'lorentzianPsi',0.5,100);
resultTuk100 = aniso(im,'tukeyPsi',0.5,100);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor100,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk100,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin100(:,32) resultLor100(:,32) resultTuk100(:,32)]);

% Lots o' iterations.  Tukey converges wheres lorentzian and
% linear continue to blur more and more.
resultLin1000 = aniso(im,'linearPsi',0,1000);
resultLor1000 = aniso(im,'lorentzianPsi',0.5,1000);
resultTuk1000 = aniso(im,'tukeyPsi',0.5,1000);
figure(1); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultLor1000,[0,3]);
figure(2); clf; 
subplot(1,2,1); imshow(im,[0,3]); subplot(1,2,2); imshow(resultTuk1000,[0,3]);
figure(3); clf; 
plot([im(:,32) resultLin1000(:,32) resultLor1000(:,32) resultTuk1000(:,32)]);

% Tradeoff between lambda and iterations.  All of these should
% give the same result
result5 = aniso(im,'tukeyPsi',0.5,5,0.5);
result10 = aniso(im,'tukeyPsi',0.5,10,0.25);
result100 = aniso(im,'tukeyPsi',0.5,100,0.025);
result1000 = aniso(im,'tukeyPsi',0.5,1000,0.0025);
figure(1); clf;
plot([im(:,32) result5(:,32) result10(:,32) result100(:,32) result1000(:,32)]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,x,Icor,B0] = BiasCorrLEMSS2D(I,Imask,V_mean,V_std,options,B0);
% BiasCorrLEMSS2D: Correct for Inhomo by local entropy minimization with splines support
%
%   [B,x,B0] = BiasCorrLEMSS2D(I,Imask,V_mean,V_std,options,B0);
%   B: Final bias filed
%   x: corrected iamge (x=I/B)
%   B0: Initial bias field in case it has been modified by the funciton
%   I: input image (should be in  [0 255])
%   mask: idientify the pixels subject to bias field
%   V_mean: mean of the background 
%   V_std: std of the background
%   options: structure with options. to innitialize the stucture call the
%   fucntion without input argument
% option = BiasCorrLEMSS2D
% option = 
%            Nknots: [30 30] knot spacings in y and x direction
%          NiterMax: 3 number of iteration max
%      flag_display: 1 if 1, display intermediate results in a new figure
%             overs: -1, if 1 pad the image to get an even number of knots
%         normalize: 1, if 1 notrmalize the image
%     flag_allknots: 1, if 1 optimze the knots at the image border
%        GainSmooth: 0, gain to force the spline to be smooth by constraining the second derivatives 
%             Bgain: 0.5000, average of the bias field over the pixels with signal            
%   B0: Initial Bias field
%
%   Olivier Salvado, 20-jan-04, Case Western Resrve University
%  Salvado et al. IEEE TMI 25(5):539-552


conds = 'v';

% --- define the default options
optionsdef.Nknots = [30 30];    % how meny pixels between knots
optionsdef.NiterMax = 3;    % 3 iteration max
% optionsdef.flag_filter = 0;     % no AD filter prior
optionsdef.flag_display = 1;    % all display on
optionsdef.overs = -1;       % automatic with Nknots to get even number of knots
optionsdef.normalize = 1;   % do normalize images
optionsdef.flag_allknots = 1;   % update all the knots
optionsdef.GainSmooth = 0;
% optionsdef.PolyOrder = 3;
optionsdef.Bgain = 0.5;

if nargin == 0, % no vargin, return the options
    B = optionsdef;
    x = [];
    return
end

if ~exist('options','var'), % if no option specify, use the default
    options = optionsdef;
end


% --- read the options
Nknots = options.Nknots;
NiterMax = options.NiterMax;
% flag_filter = options.flag_filter;
flag_display = options.flag_display;
overs = options.overs;
normalize = options.normalize;
flag_allknots =options.flag_allknots;
GainSmooth = options.GainSmooth;
% PolyOrder = options.PolyOrder;
Bgain = options.Bgain;

randn('state',0);
%disp('')
%disp(' -------------------- BIAS correction ------------------------')
%disp('    2D With Local Entropy Minimization with Spline support  V3')

% -- log transform
logtransform = 0;
if logtransform==0,
%    disp('  Using multiplicative bias: y=bx')
elseif logtransform==1
%    disp('  Using Log transformed data: y=x+b');
end


%clf % zeynep
%set(gcf,'BackingStore','on');

% --- add 'overs' lines with rt  around
[nr,nc] = size(I);
if overs<0, % automatic wiht Nknots,
    temp = nr/Nknots(1);
    k = floor(temp);
    epsilon = temp-k;
    nrp = Nknots(1)*(k+1);  % new dimension
    oversr = ceil((nrp-nr)/2);   % overs to get the new diimension
    
    temp = nc/Nknots(2);
    k = floor(temp);
    epsilon = temp-k;
    ncp = Nknots(2)*(k+1);  % new dimension
    oversc = ceil((ncp-nc)/2);   % overs to get the new diimension
    
else
    oversr = overs;
    oversc = overs;
end

if overs~=0, % need to oversize
    Imask = padarray(Imask,[oversr oversc],0,'both');
    I = padarray(I,[oversr oversc],V_mean,'both');
    if exist('B0','var'),
        B0 = padarray(B0,[oversr oversc],'replicate');
    end
end

% --- normalize the image
if normalize, % no normalization to be able to compare B to Btrue...
    If = imfilter(I,fspecial('gaussian',31,11));
    Rnorm = 100/max(If(:));
    y = double(I)*Rnorm;
    V_mean = V_mean * Rnorm;
    V_std = V_std * Rnorm;
else
    y = double(I);
    Rnorm = 1;
end

% --- Initialize B
mask = (Imask == 1);
Gauss = fspecial('gaussian',30,10);
pol = 0;
pol_old = 1;

% --- Initialize B
mask = (Imask == 1);
Gauss = fspecial('gaussian',30,10);
pol = 0;
pol_old = 1;

% --- build the grid
i = [1:size(y,1)];
j = [1:size(y,2)];
[ii,jj] = ndgrid(i,j);

% --- get the knots and B(knots)
NknotsI = round(size(y,1)/Nknots(1));
NknotsJ = round(size(y,2)/Nknots(2));
ik = round(linspace(1,size(y,1),NknotsI));
jk = round(linspace(1,size(y,2),NknotsJ));
[iik,jjk] = ndgrid(ik,jk);

% --- tag the knots to be updated
if 0,
    Mknots = zeros(NknotsI,NknotsJ);
    for ki=1:NknotsI,
        for kj=1:NknotsJ,
            if mask(ik(ki),jk(kj)) == 1,
                Mknots(ki,kj) = 1;
            end
        end
    end
    Mknots = (Mknots>0);
else,   % make all the knots valid except the borders
    Mknots = ones(NknotsI,NknotsJ);
    Mknots(1,:) = 0;
    Mknots(end,:) = 0;
    Mknots(:,1) = 0;
    Mknots(:,end) = 0;
end


% ik,jk: all the knots vectors size(ik)=[1,Nknots];
% Mknots: 1 when the knots need to be optimized
Bk = B0(ik,jk);         % B at the knots

% --- sort the knots from higher B to lower B
[dummy,knotlist] = sort(-Bk(:).*Mknots(:));       % go from brightest to dimmest
% knotlist = [1:Nknots*Nknots];                                   % does not sort the knots
Bk = Bk + 0*randn(size(Bk));
Bk0 = Bk;
Hx = inf;

% --- compute for B0
pp0 = csape({ik,jk}, Bk0,conds);
B0 = fnval(pp0,{i,j});
% ratioB0 = max(B0(mask));
ratioB0 = mean(B0(mask))/Bgain;
B0 = B0/ratioB0;
x0 = y;
x0 = y./B0.*Imask + (1-Imask).*y;

% --- main loop
Nktot = NknotsI*NknotsJ;
tic
iter = 1;
again = 1;
B = B0;
B = Bgain*B/mean(B(mask));

while again;   % iteration of the fitting
    Mopt = zeros(size(y));
    Mlead = zeros(size(y));
    Bkold = Bk;
    Hxold = Hx;
    for idx=1:Nktot,           % iteration on the knots
        k = knotlist(idx);        % next on the list
        %         k = idx;                        % does not sort the knots
        if (flag_allknots | (Mknots(k)>0)),          % this knot needs to be updated (all the knots now)
            
            if flag_display,
                subplot(243)
                hold on, plot(jjk(k),iik(k),'r+'), hold off
                subplot(244)
                hold on, plot(jjk(k),iik(k),'r+'),
                plot(jjk(:),iik(:),'co')
                hold off
                drawnow
            end
            
            % --- get the incremental mask around the knot
            Mnew = zeros(size(y));          % blank mask
            if (iik(max(k-1,1))<iik(k)), imin = iik(k-1); else imin=iik(k);end
            if (iik(min(k+1,Nktot))>iik(k)), imax = iik(k+1); else imax=iik(k);end
            if (jjk(max(k-NknotsI,1))<jjk(k)), jmin = jjk(k-NknotsI); else jmin=jjk(k);end
            if (jjk(min(k+NknotsI,Nktot))>jjk(k)), jmax = jjk(k+NknotsI); else jmax=jjk(k);end
            zonei = (ii>=imin) & (ii<=imax);              % mask for all the i
            zonej = (jj>=jmin) & (jj<=jmax);    % mask for all the j
            Mnew = zonei & zonej & mask;
            
            % --- get the area within 2 knots around the knot under
            % optimization
            zone = 5;
            whik = mod( k , length(ik) );   
            whjk = ceil( k / length(ik) );
            imin = ik( max( 1 ,whik-zone) );
            imax = ik( min( NknotsI ,whik+zone) );
            jmin = jk( max( 1 ,whjk-zone) );
            jmax = jk( min( NknotsJ ,whjk+zone) );
            area =[imin imax jmin jmax];    
            
            if sum(Mnew(:))>300, % There should be enough data in the mask
                
                if idx<6000,
                    Mlead = Mlead | Mnew;
                end
                Mopt = Mlead | Mnew;
                
                % --- optimize entropy by moving the knot k
                % ============================
                option = optimset('TolX',1e0,'Display','off',...
                    'DiffMaxChange',1,...
                    'DiffMinChange',0.01,...
                    'NonlEqnAlgorithm','lm');
                %                 'DiffMaxChange',max(6-idx,1),...
                %                 'DiffMinChange',1/(idx+1)^2,...
                
                %     Bkoptim = fminsearch(@optfungo0,Bk(k),option,ik,y,Bk,k,i);
                upperBnd = Bk(k)*1.15;
                lowerBnd = max( 0.1 , Bk(k)*0.85 );
                [Bkoptim] = fminbnd(@optfungo7,lowerBnd,upperBnd,option,...
                    ik,jk,y,Bk,k,i,j,Mopt,mask,GainSmooth,area,B,Bgain);
                
                if Bkoptim>0,
                    Bk(k) = Bkoptim;
                    
                    % --- find the spline
                    %     pp = csapi(ik, Bk, i);
                    pp = csape({ik,jk}, Bk,conds);
                    
                    % --- interpolate B
                    %     B = fnval(spi,i);
                    B = fnval(pp,{i,j});
                    %                     ratioB = max(B(mask));
                    ratioB = mean(B(mask))/Bgain;
                    B = B/ratioB;
                    
                    % --- reconstruct 
                    x = y;
                    %                     x(mask) = y(mask)./B(mask);
                    x = y./B.*Imask + (1-Imask).*y;

                end
                PDFx = hist(x(mask),[1:1:300]);                
                % --- display
                if flag_display,
                    subplot(221)
                    row = round(nr/2);
%                     ycor = y./B;
%                     ycor(~mask) = y(~mask);
%                     ycor = y./B.*Imask + (1-Imask).*y;
                    
                    plot([B(row,:)*200 ; B0(row,:)*200 ; y(row,:)  ; x(row,:)]')
                    legend('B estimated','B0','y','xhat')
                    ylim([0 300])
                    
                    subplot(243)
                    imagesc(x),axis image,title('image reconstructed x'),axis off
                    
                    subplot(244)
                    temp = B;
                    imagesc(B,[0 1.3]),axis image,axis off,title('Bias estimated')
                    
                    subplot(248), 
%                     imagesc(Bk),axis image,axis off
                    imagesc(Mopt+Mnew),axis image, axis off
                    
                    subplot(247),imagesc(y),axis image,title('image initial'),axis off
                    %                 subplot(248),imagesc(B./Btrue,[0.9 1.1]),axis image,axis off,colorbar
                    %                 title('B/Btrue')
                    %                subplot(248),imagesc(Mopt),axis image, axis off
                    
                    subplot(223)
                    subplot(223),plot([1:1:300],PDFx), title('histogram of x estimated'),
                end                    
                    % --- entropy calculation
                    PDFx = PDFx((PDFx>0));
                    PDFx = PDFx / sum(PDFx(:));
                    Hx = -sum( PDFx.*log(PDFx) );
            end
        end % if
    end % for
    if flag_display & 1,
        disp([' Entropy: ' num2str(Hx)])
    end
    Bchange = mean( abs(Bk(:)-Bkold(:))./abs(Bkold(:)) )*100;
    Hchange = (Hx-Hxold)./abs(Hxold(:)) *100;
    again = (iter<=1) | ( (iter<NiterMax)  & (Bchange>0.5) & (abs(Hchange)>0.001));
    disp([' Change on knots:' num2str(Bchange) ' %'])
    disp([' Change on H:' num2str(Hchange) ' %'])
    iter = iter +1;
end


%ProcTime = toc

pause(1)
%clf % zeynep
if flag_display,
    ymax = max(y(:));
    subplot(221),imagesc(y,[0 ymax]),axis image,axis off,title('Original')
    subplot(222),imagesc(x,[0 ymax]),axis image,axis off,title('Corrected with LEMS')
    subplot(223),hist(y(:),[0:300]);
    subplot(224),hist(x(:),[0:300]);
    colormap(gray(256))
end

x = x / Rnorm;

% --- correct original (without filtering) for B
Icor = I./B.*Imask + (1-Imask).*I;

% --- resize to original size
if overs~=0, % need to resize
    x = x(oversr+1:end-oversr,oversc+1:end-oversc);
    Icor = Icor(oversr+1:end-oversr,oversc+1:end-oversc);
    B = B(oversr+1:end-oversr,oversc+1:end-oversc);
    B0 = B0(oversr+1:end-oversr,oversc+1:end-oversc);
end

% ------------------------ FUNCTION FOR OPTIM
function [Hx,B,x,PDFx]= optfungo7(Bkoptim,ik,jk,y,Bk,k,i,j,Mopt,M,GainSmooth,area,B,Bgain);
conds = 'v';
% --- find the spline
%     pp = csapi(ik, Bk, i);
Bk(k) = Bkoptim;
pp = csape({ik,jk}, Bk,conds);

% --- interpolate B
% B = fnval(pp,{i,j});
% % B = B/max(B(M));
% B = 0.5*B/mean(B(M));

if 0,
    B = fnval(pp,{i,j});
else
    Barea = fnval(pp,{i(area(1):area(2)),j(area(3):area(4))});
    B(area(1):area(2),area(3):area(4)) = Barea;
end
B = Bgain*B/mean(B(M));

% --- reconstruct 
x = y./B;
% x(~M) = y(~M);  % the bias does not affect the background

% --- entropy calculation
PDFx = hist(x(Mopt),[1:0.5:400]);
% PDFx = filtfilt([0.5],[1 -0.5],PDFx);
PDFx = PDFx((PDFx>0));
Hx = -sum( PDFx.*log(PDFx) );

if GainSmooth>0   % use a smoothness constraint
    N = length(y(Mopt));
    Hmax= log(N);
    Hmin = -N*log(N);
    Hx = -(Hx - Hmax)/Hmin;             % normalized entropy
    
    if 0,
        Bder2 = fnder(pp,[1 1]');
        Bder2 = (fnval(Bder2,{i,j})/N).^2;
        Smoothness = sum(Bder2(:));
    else
        Eb2 = imfilter(Bk,[0 -1 0;-1 4 -1;0 -1 0],'full','replicate');
        Eb4 = imfilter(Eb2,[0 -1 0;-1 4 -1;0 -1 0],'full','replicate');
        Eb4 = (Eb4(3:end-2,3:end-2)/N).^2;
        Smoothness = sum(Eb4(:));
    end
    Hx = Hx + GainSmooth*Smoothness;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,basis,p]=PolyMaskFilter(in,order,mask,basis)
% PolyMaskFilter:   Polynomial filter under mask constraint
%
%   [out,basis,p]=PolyMaskFilter(in,order,mask,basis)
%   Poly filter: fit a ORDERth polynomial function on the data and
%       interpolate on all the image
%
% CWRU, OS 05mar03
%  modif by OS may03 to use weighted LS: cov matrix too big for an image
%                   does not work
%  modif by OS 29may03, use a subsample if Weights are given and implement
%                   weighted least square
%  modif by OS 22sep03, use orthogonal polynomial. Can be given as an input
%                   argument,... not so ppolymial


% --- where to fit the polynomial function
if ~exist('mask'),
    mask = (abs(in)>0.0001);
else
    mask = logical(mask);
end

% --- test if weighted LS 
isw = sum(sum(mask>0 & mask<1))>1;
DIM = size(in);

% --- compute the basis
if ~exist('basis')
%     flag_basis = 0;
    Ni = DIM(1)*DIM(2);
    [x,y] = meshgrid(0:DIM(1)-1,0:DIM(2)-1);
    x = x'/(DIM(1)-1) - 0.5;
    y = y'/(DIM(2)-1) - 0.5;
    
    % --- build the basis
    kp = 0;
    for kx=0:order,
        for ky=0:order,
            if (kx+ky)<=order,
                kp = kp+1;
                basis(:,:,kp) = x.^kx .* y.^ky;
%                 [kx ky]
            end
        end
    end
end

% newI    = in.*(mask>0);
newI = in;
newI(~mask) = 0;
[wx,wy] = find( newI > 0);
idx     = sub2ind(size(in),wx,wy);
Ydata   = in(idx);
W = mask(idx);

X = [];
for k=1:size(basis,3),
    temp = basis(:,:,k);
    X = [X temp(idx)]; 
end

    
if isw,
    Ns = length(Ydata);
    % --- get a subsample of size 1000 max (about 1000)
    if Ns>1000,
        di = round(Ns/1000);
        index = [1:di:Ns]';
%         Xdata = Xdata(index,:);
        X = X(index,:);
        W = diag(W(index));
        Ydata = Ydata(index);
    end
end

if ~isw,
    p = X\Ydata;
else
    p = inv(X'*W*X)*X'*W*Ydata;
end

F = reshape(basis,DIM(1)*DIM(2),size(basis,3))*p;
F = reshape(F,DIM(1),DIM(2));
out = F;

