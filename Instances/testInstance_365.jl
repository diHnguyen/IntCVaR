edge = [1 3; 1 7; 1 13; 1 15; 2 3; 2 6; 2 7; 2 9; 2 10; 2 14; 3 2; 3 5; 3 8; 3 9; 3 14; 3 15; 4 5; 4 14; 4 15; 5 3; 5 4; 5 7; 5 8; 6 3; 6 8; 6 11; 6 12; 6 13; 6 14; 7 2; 7 5; 7 8; 7 9; 7 12; 7 13; 7 14; 8 4; 8 9; 8 10; 8 13; 8 14; 9 3; 9 7; 9 10; 9 11; 9 15; 10 2; 10 3; 10 7; 10 8; 10 9; 10 13; 11 2; 11 3; 11 9; 11 10; 11 14; 12 4; 12 9; 12 11; 12 15; 13 4; 13 6; 13 9; 13 11; 13 12; 14 4; 14 8; 14 10; 14 13; 14 15]
cL_orig = [71.0, 37.0, 124.0, 18.0, 46.0, 109.0, 203.0, 341.0, 183.0, 294.0, 17.0, 45.0, 161.0, 167.0, 136.0, 225.0, 35.0, 104.0, 28.0, 4.0, 12.0, 20.0, 98.0, 24.0, 9.0, 233.0, 33.0, 39.0, 11.0, 46.0, 73.0, 27.0, 31.0, 93.0, 104.0, 255.0, 151.0, 42.0, 64.0, 60.0, 56.0, 205.0, 80.0, 24.0, 31.0, 173.0, 289.0, 175.0, 89.0, 33.0, 4.0, 44.0, 79.0, 7.0, 5.0, 49.0, 53.0, 313.0, 90.0, 49.0, 26.0, 235.0, 39.0, 23.0, 34.0, 6.0, 175.0, 211.0, 39.0, 27.0, 6.0]
cU_orig = [71.0, 37.0, 124.0, 128.0, 46.0, 169.0, 203.0, 341.0, 273.0, 294.0, 17.0, 45.0, 161.0, 167.0, 136.0, 225.0, 35.0, 780.0, 28.0, 20.0, 12.0, 20.0, 98.0, 24.0, 9.0, 233.0, 33.0, 39.0, 11.0, 46.0, 73.0, 35.0, 159.0, 167.0, 104.0, 433.0, 151.0, 42.0, 64.0, 60.0, 220.0, 205.0, 80.0, 26.0, 31.0, 173.0, 289.0, 175.0, 89.0, 33.0, 4.0, 76.0, 79.0, 103.0, 21.0, 49.0, 53.0, 313.0, 90.0, 49.0, 26.0, 519.0, 39.0, 23.0, 34.0, 20.0, 175.0, 211.0, 145.0, 27.0, 46.0]
d = [19.0, 36.0, 7.0, 47.0, 13.0, 27.0, 49.0, 42.0, 10.0, 37.0, 12.0, 17.0, 14.0, 20.0, 3.0, 35.0, 42.0, 17.0, 1.0, 30.0, 46.0, 41.0, 10.0, 8.0, 50.0, 15.0, 3.0, 24.0, 11.0, 38.0, 7.0, 22.0, 38.0, 32.0, 15.0, 33.0, 2.0, 40.0, 43.0, 43.0, 23.0, 19.0, 10.0, 13.0, 36.0, 13.0, 45.0, 49.0, 39.0, 15.0, 40.0, 44.0, 2.0, 28.0, 34.0, 40.0, 39.0, 40.0, 21.0, 38.0, 33.0, 14.0, 24.0, 41.0, 28.0, 33.0, 8.0, 14.0, 15.0, 4.0, 20.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
