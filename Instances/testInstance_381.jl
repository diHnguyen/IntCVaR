edge = [1 3; 1 7; 1 11; 1 14; 2 3; 2 7; 2 9; 2 10; 2 11; 2 15; 3 2; 3 6; 3 8; 3 9; 3 11; 3 13; 3 14; 3 15; 4 2; 4 3; 4 8; 4 11; 4 12; 5 8; 5 9; 5 10; 5 11; 5 13; 5 15; 6 2; 6 3; 6 9; 6 10; 7 4; 7 6; 7 9; 7 13; 7 15; 8 5; 8 6; 8 10; 8 11; 8 15; 9 6; 9 7; 9 12; 9 14; 9 15; 10 6; 10 15; 11 5; 11 13; 11 14; 12 5; 12 7; 12 9; 12 10; 12 15; 13 4; 13 9; 13 11; 13 12; 14 2; 14 11; 14 12; 14 15]
cL_orig = [83.0, 138.0, 57.0, 639.0, 38.0, 123.0, 88.0, 303.0, 169.0, 33.0, 17.0, 21.0, 157.0, 28.0, 37.0, 273.0, 476.0, 561.0, 29.0, 9.0, 124.0, 27.0, 323.0, 93.0, 17.0, 120.0, 214.0, 159.0, 334.0, 122.0, 108.0, 94.0, 199.0, 86.0, 45.0, 39.0, 130.0, 350.0, 32.0, 43.0, 43.0, 60.0, 13.0, 47.0, 66.0, 30.0, 4.0, 110.0, 192.0, 203.0, 20.0, 29.0, 29.0, 14.0, 12.0, 93.0, 24.0, 0.0, 319.0, 21.0, 6.0, 22.0, 432.0, 121.0, 60.0, 5.0]
cU_orig = [83.0, 138.0, 57.0, 639.0, 54.0, 127.0, 214.0, 463.0, 169.0, 543.0, 17.0, 21.0, 203.0, 28.0, 45.0, 273.0, 476.0, 561.0, 29.0, 9.0, 258.0, 27.0, 323.0, 93.0, 17.0, 120.0, 350.0, 159.0, 334.0, 122.0, 108.0, 132.0, 199.0, 86.0, 45.0, 71.0, 130.0, 350.0, 32.0, 63.0, 111.0, 60.0, 13.0, 119.0, 66.0, 128.0, 16.0, 422.0, 192.0, 203.0, 20.0, 29.0, 37.0, 14.0, 374.0, 93.0, 50.0, 6.0, 319.0, 21.0, 22.0, 44.0, 728.0, 121.0, 60.0, 5.0]
d = [24.0, 48.0, 49.0, 2.0, 20.0, 15.0, 49.0, 50.0, 32.0, 19.0, 48.0, 36.0, 15.0, 42.0, 50.0, 15.0, 26.0, 33.0, 38.0, 3.0, 48.0, 24.0, 49.0, 25.0, 38.0, 23.0, 9.0, 33.0, 15.0, 31.0, 10.0, 27.0, 43.0, 3.0, 17.0, 26.0, 41.0, 41.0, 31.0, 5.0, 7.0, 33.0, 33.0, 25.0, 24.0, 5.0, 23.0, 25.0, 21.0, 9.0, 17.0, 43.0, 35.0, 15.0, 10.0, 1.0, 37.0, 37.0, 28.0, 7.0, 44.0, 46.0, 33.0, 37.0, 45.0, 42.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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