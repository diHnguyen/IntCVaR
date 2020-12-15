edge = [1 3; 1 5; 1 8; 1 9; 1 15; 2 4; 2 8; 2 9; 2 10; 2 12; 3 2; 3 5; 3 10; 3 13; 3 14; 3 15; 4 2; 4 5; 4 8; 4 13; 4 15; 5 7; 5 9; 5 10; 5 13; 5 14; 6 8; 6 11; 6 12; 6 13; 6 15; 7 12; 7 13; 7 15; 8 7; 8 10; 8 12; 8 13; 8 14; 9 2; 9 3; 9 4; 9 8; 9 11; 10 7; 10 8; 10 12; 11 5; 11 6; 11 10; 11 12; 11 15; 12 2; 12 8; 12 10; 13 2; 13 6; 13 9; 13 10; 13 11; 13 12; 13 14; 13 15; 14 7; 14 8; 14 9]
cL_orig = [15.0, 193.0, 50.0, 41.0, 220.0, 84.0, 50.0, 88.0, 184.0, 454.0, 6.0, 55.0, 107.0, 156.0, 526.0, 94.0, 63.0, 2.0, 91.0, 60.0, 501.0, 11.0, 78.0, 112.0, 51.0, 443.0, 13.0, 103.0, 184.0, 326.0, 99.0, 3.0, 55.0, 335.0, 16.0, 7.0, 9.0, 84.0, 289.0, 27.0, 235.0, 237.0, 1.0, 95.0, 75.0, 49.0, 83.0, 182.0, 95.0, 29.0, 39.0, 25.0, 15.0, 99.0, 59.0, 37.0, 332.0, 145.0, 8.0, 12.0, 18.0, 3.0, 84.0, 2.0, 236.0, 221.0]
cU_orig = [71.0, 193.0, 66.0, 41.0, 220.0, 84.0, 102.0, 88.0, 184.0, 454.0, 34.0, 55.0, 107.0, 156.0, 526.0, 94.0, 63.0, 28.0, 91.0, 104.0, 501.0, 11.0, 78.0, 112.0, 51.0, 443.0, 139.0, 103.0, 184.0, 326.0, 99.0, 59.0, 55.0, 335.0, 54.0, 35.0, 9.0, 110.0, 289.0, 27.0, 235.0, 237.0, 55.0, 95.0, 75.0, 49.0, 83.0, 182.0, 95.0, 29.0, 39.0, 231.0, 15.0, 99.0, 119.0, 37.0, 332.0, 145.0, 26.0, 46.0, 18.0, 19.0, 84.0, 6.0, 236.0, 249.0]
d = [17.0, 29.0, 40.0, 2.0, 11.0, 34.0, 18.0, 29.0, 43.0, 13.0, 26.0, 49.0, 21.0, 34.0, 21.0, 10.0, 10.0, 46.0, 27.0, 12.0, 27.0, 41.0, 22.0, 38.0, 26.0, 48.0, 17.0, 50.0, 39.0, 28.0, 33.0, 2.0, 34.0, 38.0, 30.0, 31.0, 6.0, 47.0, 38.0, 9.0, 25.0, 42.0, 13.0, 50.0, 11.0, 14.0, 30.0, 39.0, 18.0, 14.0, 22.0, 39.0, 11.0, 9.0, 45.0, 11.0, 2.0, 37.0, 36.0, 35.0, 2.0, 4.0, 42.0, 40.0, 35.0, 34.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
