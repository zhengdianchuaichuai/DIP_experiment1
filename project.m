%% Play video
v = VideoReader('Calcium500frames.avi');
currAxes = axes;
while hasFrame(v)
    vidFrame = readFrame(v);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v.FrameRate);
    disp(v.CurrentTime)
end
%% Load variables
load('matlab.mat');
%% Show single image
imshow(frames(:,:,50));
%% Get normal pd for histogram data
pd_normal = fitdist(reshape(frames(:,:,:),[],1),'Normal');
%% Get 
%% Hist of all images
figure;
hold on;
[counts,binlocs] = imhist(frames(:,:,:));

y_normal = pdf(pd_normal,binlocs);
y_poiss = poisspdf((binlocs),15);

bar(binlocs,counts/sum(counts));
plot(binlocs,y_normal);
plot(binlocs,y_poiss)

legend('data','normal','poisson');

hold off;
%% Get DFT of noise
sample = frames(1:250,1:250,50);
DFT_im(sample);