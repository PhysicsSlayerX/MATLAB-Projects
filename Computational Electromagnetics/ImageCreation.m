%Image Creation

close all;
clc;
clear all;


figure
image=zeros(360,360,3); %initialize
image(:,:,1)=rand(360,360);   %Red (dark red)
image(:,:,2)=rand(360,360);   %Red (maximum value)
image(:,:,3)=rand(360,360);  %Green
image(:,:,1)=rand(360,360);
imshow(image)
