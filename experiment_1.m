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
% % 
% % 
% % 
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
In this dataset: When a neuron is fully activated, the diameter could be up to ~ 30 pixels. 

So 4-30 pixels is considered as neuron activation size.

------ Logic: ----- 
1) set up intesity threshold F0 and remove noise whose intensity is lower than 1.2F0 (by changing the image into binary)
2) set up size thredshold to find out the neuron (size: 5-30 pixels due to limitation of 'imfindcircles'; shape: circle)
3) loop thru all frames in order and count firing event by detecting black -> white

%}


F0 = 55;
delta_F = 1.2*F0;
centersPrevFrame = zeros(1,2);
centersPrevPrevFrame = zeros(1,2);
num_of_circles = 0;
centersSum = {};
radiiSum = {};
is = [];


for i = 1:500
    fileName = strcat('frame', num2str(i), '.jpg');
    img = rgb2gray(imread(sprintf('frame%d.jpg',i)));
%     se = strel('disk',1);
%     imgDenoised = imerode(img,se); % erode
    

%%% step 1) %%%
    img(img<delta_F) = 0;
    img(img>=delta_F) = 255;
%     figure; imagesc(img); colormap gray;

%%% step 2) & 3) %%%
    [centers,radii] = imfindcircles(img,[5 30],'Sensitivity',0.86);
    if ~isempty(centers)
%         figure; imshow(img);
%         h = viscircles(centers,radii);
        for n = 1:size(centers,1)
            dist1 = sqrt((centers(n,1) - centersPrevFrame(:,1)).^2 + ((centers(n,2) - centersPrevFrame(:,2)).^2)); % distance to centers in the previous frame
            dist2 = sqrt((centers(n,1) - centersPrevPrevFrame(:,1)).^2 + ((centers(n,2) - centersPrevPrevFrame(:,2)).^2)); % distance to centers in the previous previous frame
            flag = 0;
            for p = 1:size(dist1)
                if find(abs(dist1(p)) < 2)
                    flag = 1;
                    break
                end
            end
            if flag == 0
                for q = 1:size(dist2)
                    if find(abs(dist2(q)) < 2)
                        break
                    end
                end
                num_of_circles = num_of_circles + 1;
                centersSum{end+1} = centers(n,:);
                radiiSum{end+1} = radii(n,:);
                is = [is i];
            end  
        end
        centersPrevPrevFrame = centersPrevFrame;
        centersPrevFrame = centers;
    end
% %     if ismember(i,[20 30 40 50 60])
% %         if isequal(centers, [0,0])
% %             figure; imshow(img); title(['Frame',num2str(i)]);
% %         else
% %             figure; imshow(img); title(['Frame',num2str(i)]);
% %             h = viscircles(centers,radii);
% %         end
% %     end
end

%% Q4

value = 0;
values = [];
for j = 1:length(is)
    i = is(j);
    img = rgb2gray(imread(sprintf('frame%d.jpg',i)));
    value = double(img(uint16(centersSum{j}(2)),uint16(centersSum{j}(1))))...
        + double(img(uint16(centersSum{j}(2)) - 1,uint16(centersSum{j}(1)) - 1))...
        + double(img(uint16(centersSum{j}(2)) + 1,uint16(centersSum{j}(1)) - 1))...
        + double(img(uint16(centersSum{j}(2)) - 1,uint16(centersSum{j}(1)) + 1))...
        + double(img(uint16(centersSum{j}(2)) + 1,uint16(centersSum{j}(1)) + 1));
    value = value/5.0;
    value = (value - F0)/F0;
    values = [values value];
end
values(values<0) = 0;
histogram(values,10);
xlabel("delF/F0");
ylabel("counts");
title("Histogram of relative intensities over firing events");

%% Q5

imgSum = zeros(512,512);
figure; imshow(imgSum);
hold on

% plot circles
for k = 1:num_of_circles
    centerK = centersSum{k};
    radiiK = radiiSum{k};
    theta = 0 : (2 * pi / 10000) : (2 * pi);
    pline_x = radiiK * cos(theta) + centerK(1);
    pline_y = radiiK * sin(theta) + centerK(2);
    plot(pline_x, pline_y, '-');
end


% count circles and display numbers
countSum{1,1} = centersSum{1}; % countSum - a cell to record the centers and count for all firing events
countSum{1,2} = 1;
for m = 2:num_of_circles
    flagS = 0;
    for s = 1:size(countSum,1)
        dist = sqrt((centersSum{m}(1) - countSum{s}(1)).^2 + (centersSum{m}(2) - countSum{s}(2)).^2);
        if dist <= 5
            countSum{s,2} = countSum{s,2} + 1;
            flagS = 1;
            break
        end
    end
    
    if flagS == 0
        countSum{s+1,1} = centersSum{m};
        countSum{s+1,2} = 1;
    end
end
        

for t = 1:size(countSum,1)
    num = num2str(countSum{t,2});
    str = strcat('\leftarrow',num);
    text(countSum{t,1}(1)+8,countSum{t,1}(2),str,'Color','red','FontSize',14)
end
hold off








