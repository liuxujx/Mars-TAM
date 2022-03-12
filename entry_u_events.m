function [value,isterminal,direction]=entry_u_events(t,inp) %
global Vf

% r  = inp(1);
V  = inp(2);
% gamma  = inp(3);
value = V-Vf;
isterminal=1;
direction=0; % -1 value´ÓÕýµ½0
end