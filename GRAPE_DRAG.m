t=-5:0.1:5;
sm=1;
tg=2;

y1=1.5232*exp(-(t-tg/2).^2/2/sm^2);
y2=tanh(t/sm)-tanh((t-tg)/sm);

plot(t,y1,t,y2)