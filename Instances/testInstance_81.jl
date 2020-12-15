edge = [1 2; 1 3; 1 4; 1 7; 1 9; 1 11; 1 12; 1 15; 2 5; 2 6; 2 9; 2 10; 2 14; 3 6; 4 5; 4 10; 4 13; 5 3; 5 9; 5 10; 5 11; 5 13; 5 15; 6 2; 6 5; 6 8; 6 9; 7 2; 7 14; 8 2; 8 9; 8 12; 8 13; 8 14; 8 15; 9 3; 9 6; 9 12; 9 14; 10 2; 10 4; 10 5; 10 6; 10 8; 10 11; 10 12; 11 7; 11 15; 12 5; 12 6; 12 9; 12 10; 12 13; 12 14; 13 2; 13 14; 14 4; 14 5; 14 8; 14 10; 14 15]
cL_orig = [34.0, 99.0, 62.0, 144.0, 206.0, 200.0, 102.0, 271.0, 150.0, 11.0, 227.0, 24.0, 196.0, 36.0, 8.0, 185.0, 91.0, 85.0, 159.0, 165.0, 197.0, 109.0, 227.0, 200.0, 10.0, 68.0, 11.0, 71.0, 341.0, 50.0, 2.0, 97.0, 127.0, 139.0, 124.0, 207.0, 40.0, 51.0, 153.0, 186.0, 32.0, 92.0, 151.0, 1.0, 19.0, 5.0, 21.0, 22.0, 284.0, 252.0, 7.0, 25.0, 27.0, 59.0, 482.0, 1.0, 82.0, 168.0, 58.0, 140.0, 38.0]
cU_orig = [34.0, 99.0, 148.0, 144.0, 206.0, 200.0, 192.0, 271.0, 150.0, 79.0, 227.0, 86.0, 196.0, 36.0, 8.0, 185.0, 91.0, 85.0, 159.0, 165.0, 197.0, 109.0, 227.0, 200.0, 10.0, 68.0, 19.0, 191.0, 341.0, 50.0, 14.0, 97.0, 127.0, 301.0, 124.0, 207.0, 40.0, 51.0, 153.0, 186.0, 168.0, 92.0, 151.0, 23.0, 19.0, 5.0, 189.0, 22.0, 284.0, 252.0, 203.0, 115.0, 27.0, 59.0, 492.0, 55.0, 82.0, 168.0, 58.0, 140.0, 38.0]
d = [38.0, 12.0, 33.0, 26.0, 15.0, 13.0, 43.0, 22.0, 32.0, 13.0, 17.0, 21.0, 8.0, 5.0, 35.0, 45.0, 47.0, 39.0, 11.0, 24.0, 25.0, 7.0, 46.0, 34.0, 9.0, 14.0, 18.0, 12.0, 43.0, 47.0, 43.0, 35.0, 21.0, 2.0, 18.0, 38.0, 11.0, 21.0, 20.0, 30.0, 2.0, 40.0, 36.0, 44.0, 1.0, 17.0, 31.0, 34.0, 26.0, 46.0, 15.0, 46.0, 8.0, 14.0, 4.0, 1.0, 19.0, 13.0, 5.0, 1.0, 30.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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