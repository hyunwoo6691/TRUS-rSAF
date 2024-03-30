clc; clear all; close all;

%%
radius = 5e-3;
height = 5e-3;

num_syn = 4;
delta_alpha = (60/127)*pi/180;

side_angle = atand(0.5*height/radius)*pi/180;

total_angle = delta_alpha*(num_syn-1)+2*side_angle;

aperture_size = radius * total_angle 

%%

num_syn = 32;
delta_alpha = (60/127)*pi/180;
sweep_angle = delta_alpha*(num_syn-1)*180/pi; % degree

r_p = sqrt(radius^2 + (0.5*height)^2);
angle_ = sweep_angle + atand(0.5*height/radius);

aperture_size = 2*r_p*sind(angle_)