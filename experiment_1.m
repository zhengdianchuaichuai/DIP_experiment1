close all;clear;clc

% read the .avi, extract the frames

% % mov = VideoReader('Calcium500frames.avi');
% % for frameNum = 1:mov.NumberOfFrames
% %     filename = strcat('frame',num2str(frameNum),'.jpg');
% %     img = read(mov, frameNum);
% %     imshow(img);
% %     imwrite(img,filename);
% % end



%% Q1
%{
These images are noisy. Characterize the noise in the video sequence. 
Feel free to use information from the 2-photon microscopy literature 
in addition to extracting statistics from the actual images.
%}


% % img1 = imread('frame1.jpg');
% % img1 = rgb2gray(img1);
% % figure; imshow(img1);
% % figure; img1Hist = histogram(img1);




%% Q2 - denoising method #1: Wavelet filtering; denoising method #2: Wiener filter
%{
Come up with an idea to denoise these images. 
Show results for frames 50 and 100 in your report (where the first frame is frame 1). 
Compare to another denoising method. Use a measure of success to show which method is superior.
%}

% % img50 = imread('frame50.jpg');
% % img50 = rgb2gray(img50);
% % figure; imshow(img50);
% % 
% % % % % % % not good performance % % % % %
% % % % % img50Denoise1 = wdenoise2(img50, 9); % Wavelet filter 
% % % % % SNR1Frame50 = 20*log10(norm(img50(:))/norm(img50(:)-img50Denoise1(:)));
% % 
% % se = strel('disk',1);
% % img50Denoise3 = imerode(img50,se); % erode
% % img50Denoise2 = wiener2(img50,[5 5]); % Wiener filter
% % SNR3Frame50 = 20*log10(norm(double(img50(:)))/norm(double(img50(:))-double(img50Denoise3(:)))); % quantify Erode filter difference from original img
% % SNR2Frame50 = 20*log10(norm(double(img50(:)))/norm(double(img50(:))-double(img50Denoise2(:)))); % quantify Wiener filter difference from original img
% % 
% % figure;
% % subplot(1,3,1);imshow(img50);title('Noisy Frame 50');
% % subplot(1,3,2);imshow(img50Denoise3);title('Erosion Denoised Frame 50');
% % subplot(1,3,3);imshow(img50Denoise2);title('Wiener Filter Denoised Frame 50');




% % 
% % img100 = imread('frame100.jpg');
% % img100 = rgb2gray(img100);
% % % % % img100Denoise1 = wdenoise2(img100, 9); % Wavelet filter
% % % % % SNR1Frame100 = 20*log10(norm(img100(:))/norm(img100(:)-img100Denoise1(:)));
% % se = strel('disk',1);
% % img100Denoise3 = imerode(img100,se); % erode
% % img100Denoise2 = wiener2(img100,[5 5]); % Wiener filter
% % SNR3Frame100 = 20*log10(norm(double(img100(:)))/norm(double(img100(:))-double(img100Denoise3(:))));
% % SNR2Frame100 = 20*log10(norm(double(img100(:)))/norm(double(img100(:))-double(img100Denoise2(:))));
% % 
% % figure;
% % subplot(1,3,1);imshow(img100);title('Noisy Frame 100');
% % subplot(1,3,2);imshow(img100Denoise3);title('Erosion Denoised Frame 100');
% % subplot(1,3,3);imshow(img100Denoise2);title('Wiener Filter Denoised Frame 100');



%% Q3 -
%{
Come up with an algorithm to detect a neuron firing event. 
Biologists often determine a baseline level of fluorescence called F0 and then detecting spikes in intensity over baseline of 20% or more. 
This spike is called delta_F and we often measure delta_F/F0. 
Report how many firing events occurred in the movie. 
Also, show detections (if any) in frames 20, 30, 40, 50 and 60.

pixel size: 0.9 µm × 0.9 µm 
The diameter of soma of a neuron range from 3 µm to 18 µm.
In this dataset: When a neuron is fully activated, the diameter could be up
to ~ 30 pixels. 

So 4-30 pixels is considered as neuron activation size.

------ Logic: ----- 
1) set up intesity threshold F0 and remove noise whose intensity is lower than 1.2F0 (by changing the image into binary)
2) set up size thredshold to find out the neuron (size: 5-30 pixels due to limitation of 'imfindcircles'; shape: circle)
3) loop thru all frames in order and count firing event by detecting black -> white

%}


F0 = 50;
delta_F = 1.2*F0;
centersPrevFrame = [0,0];
num_of_circles = 0;


for i = 1:500
    fileName = strcat('frame', num2str(i), '.jpg');
    img = rgb2gray(imread(sprintf('frame%d.jpg',i)));

%%% step 1) %%%
    img(img<delta_F) = 0;
    img(img>=delta_F) = 255;

%%% step 2) & 3) %%%
    [centers,radii] = imfindcircles(img,[5 30],'Sensitivity',0.8);
    if isempty(centers)
        centers = [0,0];
    else
% %         figure; imshow(img);
% %         h = viscircles(centers,radii);
        dist = sqrt((centers(1,1) - centersPrevFrame(1,1))^2 + ((centers(1,2) - centersPrevFrame(1,2))^2));
        if find(abs(dist) > 2)
            num_of_circles = num_of_circles + size(centers,1);
            centersPrevFrame = centers;
        end
    end
    if ismember(i,[20 30 40 50 60])
        if isequal(centers, [0,0])
            figure; imshow(img); title(['Frame',num2str(i)]);
        else
            figure; imshow(img); title(['Frame',num2str(i)]);
            h = viscircles(centers,radii);
        end
    end
end








