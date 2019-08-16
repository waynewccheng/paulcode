data = xlsread('380to390.xlsx','B3:OL13')
plot(380:780,data')
axis([380 780 0 0.1])