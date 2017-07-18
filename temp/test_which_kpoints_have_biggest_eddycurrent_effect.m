
delta_x_shift=circshift(delta_x,1,2);
delta_y_shift=circshift(delta_y,1,2);


ints=abs(delta_y_shift)>13;
figure(81);
plot(p(1,ints,1),p(2,ints,1),'.')