edge = [1 4; 1 7; 1 8; 1 11; 1 14; 2 5; 2 11; 2 14; 2 15; 3 4; 3 5; 3 7; 3 9; 3 13; 3 14; 4 2; 4 3; 4 5; 4 8; 4 10; 4 11; 4 12; 4 14; 5 2; 5 10; 5 11; 5 15; 6 2; 6 8; 6 12; 6 13; 6 15; 7 5; 7 8; 7 9; 7 12; 8 2; 8 3; 8 4; 8 10; 8 11; 8 14; 9 3; 9 6; 9 10; 9 13; 10 2; 10 6; 10 7; 10 8; 10 12; 10 13; 10 14; 11 3; 11 6; 11 7; 11 9; 12 4; 12 8; 12 9; 12 11; 12 13; 13 2; 13 6; 13 9; 13 10; 13 12; 14 4; 14 5; 14 7; 14 12; 14 13; 14 15]
cL_orig = [120.0, 41.0, 30.0, 1.0, 91.0, 4.0, 421.0, 400.0, 172.0, 40.0, 23.0, 76.0, 223.0, 34.0, 133.0, 51.0, 31.0, 23.0, 157.0, 297.0, 199.0, 25.0, 388.0, 86.0, 44.0, 222.0, 329.0, 18.0, 65.0, 217.0, 270.0, 164.0, 62.0, 44.0, 75.0, 99.0, 29.0, 80.0, 82.0, 54.0, 54.0, 146.0, 128.0, 109.0, 45.0, 58.0, 217.0, 199.0, 45.0, 18.0, 48.0, 73.0, 133.0, 388.0, 214.0, 74.0, 10.0, 20.0, 11.0, 55.0, 34.0, 22.0, 91.0, 228.0, 67.0, 69.0, 50.0, 31.0, 2.0, 306.0, 72.0, 37.0, 2.0]
cU_orig = [120.0, 43.0, 30.0, 111.0, 91.0, 64.0, 421.0, 400.0, 1122.0, 40.0, 23.0, 76.0, 223.0, 276.0, 823.0, 51.0, 31.0, 23.0, 157.0, 297.0, 199.0, 127.0, 388.0, 86.0, 44.0, 222.0, 329.0, 90.0, 65.0, 217.0, 270.0, 164.0, 62.0, 44.0, 75.0, 351.0, 29.0, 80.0, 196.0, 54.0, 54.0, 146.0, 130.0, 109.0, 45.0, 58.0, 217.0, 199.0, 45.0, 46.0, 48.0, 73.0, 133.0, 388.0, 214.0, 74.0, 10.0, 112.0, 55.0, 55.0, 44.0, 46.0, 91.0, 228.0, 315.0, 215.0, 50.0, 31.0, 4.0, 306.0, 72.0, 37.0, 2.0]
d = [25.0, 37.0, 10.0, 5.0, 48.0, 30.0, 26.0, 15.0, 37.0, 6.0, 37.0, 33.0, 12.0, 49.0, 31.0, 42.0, 24.0, 41.0, 18.0, 30.0, 25.0, 31.0, 49.0, 29.0, 11.0, 43.0, 10.0, 1.0, 16.0, 5.0, 35.0, 19.0, 19.0, 15.0, 26.0, 15.0, 14.0, 3.0, 36.0, 24.0, 14.0, 5.0, 48.0, 24.0, 12.0, 27.0, 20.0, 10.0, 48.0, 46.0, 19.0, 2.0, 30.0, 45.0, 27.0, 7.0, 32.0, 40.0, 50.0, 41.0, 23.0, 31.0, 31.0, 12.0, 45.0, 8.0, 12.0, 36.0, 32.0, 48.0, 34.0, 14.0, 16.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
