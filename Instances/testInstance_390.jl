edge = [1 2; 1 9; 1 11; 1 12; 1 14; 2 7; 2 9; 2 10; 2 12; 2 15; 3 9; 3 13; 3 14; 3 15; 4 3; 4 6; 4 8; 4 10; 5 2; 5 13; 5 14; 6 8; 6 10; 6 12; 6 15; 7 2; 7 8; 7 10; 7 11; 7 13; 7 15; 8 3; 8 6; 8 7; 8 13; 8 14; 9 2; 9 5; 9 8; 10 2; 10 4; 10 5; 10 8; 10 12; 10 13; 10 15; 11 7; 11 8; 12 5; 12 6; 12 7; 12 10; 13 3; 13 5; 13 6; 13 8; 13 9; 13 11; 13 12; 13 14; 13 15; 14 3; 14 5; 14 9; 14 10]
cL_orig = [15.0, 80.0, 468.0, 18.0, 45.0, 89.0, 234.0, 220.0, 172.0, 204.0, 5.0, 77.0, 13.0, 204.0, 45.0, 3.0, 152.0, 164.0, 55.0, 0.0, 381.0, 86.0, 21.0, 11.0, 215.0, 87.0, 10.0, 41.0, 147.0, 126.0, 208.0, 81.0, 52.0, 1.0, 7.0, 26.0, 91.0, 81.0, 35.0, 40.0, 64.0, 92.0, 6.0, 58.0, 53.0, 42.0, 54.0, 27.0, 20.0, 3.0, 10.0, 22.0, 189.0, 199.0, 1.0, 192.0, 191.0, 46.0, 23.0, 28.0, 2.0, 337.0, 324.0, 42.0, 51.0]
cU_orig = [15.0, 176.0, 468.0, 18.0, 45.0, 89.0, 234.0, 220.0, 172.0, 204.0, 5.0, 235.0, 13.0, 340.0, 51.0, 103.0, 152.0, 164.0, 55.0, 584.0, 381.0, 86.0, 21.0, 11.0, 215.0, 87.0, 10.0, 77.0, 147.0, 126.0, 208.0, 81.0, 52.0, 1.0, 141.0, 26.0, 91.0, 81.0, 35.0, 266.0, 64.0, 92.0, 6.0, 58.0, 53.0, 278.0, 54.0, 27.0, 20.0, 13.0, 10.0, 22.0, 189.0, 199.0, 1.0, 192.0, 191.0, 46.0, 73.0, 56.0, 2.0, 337.0, 324.0, 456.0, 51.0]
d = [35.0, 7.0, 35.0, 43.0, 4.0, 16.0, 6.0, 15.0, 33.0, 29.0, 32.0, 2.0, 38.0, 30.0, 14.0, 11.0, 40.0, 16.0, 18.0, 7.0, 36.0, 5.0, 21.0, 7.0, 49.0, 8.0, 13.0, 18.0, 7.0, 47.0, 29.0, 33.0, 10.0, 27.0, 5.0, 16.0, 37.0, 16.0, 4.0, 7.0, 25.0, 45.0, 41.0, 18.0, 39.0, 40.0, 22.0, 15.0, 5.0, 38.0, 26.0, 20.0, 45.0, 7.0, 38.0, 30.0, 41.0, 10.0, 22.0, 18.0, 29.0, 38.0, 35.0, 24.0, 9.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =15

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 5
last_node = maximum(edge)
