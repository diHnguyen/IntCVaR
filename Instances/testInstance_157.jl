edge = [1 4; 1 7; 1 10; 1 11; 1 12; 2 4; 2 6; 3 2; 3 5; 3 7; 3 10; 3 11; 4 6; 4 10; 4 15; 5 3; 5 6; 5 9; 5 10; 5 11; 5 12; 6 2; 6 3; 6 8; 6 11; 7 5; 7 8; 7 11; 7 14; 8 3; 8 5; 8 6; 8 14; 8 15; 9 3; 9 4; 9 7; 10 5; 10 6; 10 8; 10 14; 11 7; 11 12; 11 13; 12 2; 12 3; 12 5; 12 6; 12 7; 12 8; 12 10; 12 11; 13 5; 13 9; 13 10; 13 14; 13 15; 14 10; 14 11]
cL_orig = [121.0, 9.0, 2.0, 231.0, 13.0, 86.0, 92.0, 48.0, 77.0, 115.0, 297.0, 340.0, 79.0, 151.0, 63.0, 59.0, 6.0, 20.0, 207.0, 160.0, 14.0, 74.0, 34.0, 34.0, 185.0, 39.0, 12.0, 199.0, 244.0, 3.0, 71.0, 34.0, 113.0, 230.0, 163.0, 216.0, 80.0, 104.0, 27.0, 51.0, 1.0, 22.0, 4.0, 2.0, 98.0, 86.0, 271.0, 85.0, 139.0, 106.0, 26.0, 21.0, 74.0, 29.0, 121.0, 2.0, 13.0, 132.0, 0.0]
cU_orig = [121.0, 561.0, 30.0, 619.0, 13.0, 102.0, 92.0, 48.0, 77.0, 115.0, 297.0, 340.0, 79.0, 151.0, 63.0, 59.0, 6.0, 20.0, 273.0, 160.0, 14.0, 74.0, 244.0, 98.0, 185.0, 39.0, 12.0, 199.0, 244.0, 263.0, 71.0, 34.0, 133.0, 230.0, 163.0, 216.0, 80.0, 104.0, 27.0, 71.0, 1.0, 96.0, 4.0, 4.0, 168.0, 102.0, 271.0, 391.0, 139.0, 106.0, 26.0, 21.0, 74.0, 29.0, 121.0, 2.0, 115.0, 132.0, 182.0]
d = [13.0, 22.0, 49.0, 25.0, 45.0, 25.0, 18.0, 28.0, 25.0, 27.0, 14.0, 13.0, 6.0, 44.0, 2.0, 33.0, 37.0, 43.0, 19.0, 17.0, 23.0, 24.0, 6.0, 7.0, 11.0, 26.0, 23.0, 6.0, 41.0, 30.0, 40.0, 34.0, 13.0, 50.0, 16.0, 11.0, 17.0, 34.0, 1.0, 4.0, 46.0, 25.0, 50.0, 37.0, 13.0, 40.0, 32.0, 31.0, 6.0, 25.0, 26.0, 31.0, 48.0, 12.0, 28.0, 45.0, 36.0, 43.0, 22.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0]

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