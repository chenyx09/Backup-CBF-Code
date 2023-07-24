function theta=round_2pi(theta)
if theta>pi
    theta=theta-2*pi*ceil(theta/2/pi);
elseif theta<-pi
    theta=theta-2*pi*floor(theta/2/pi);
end