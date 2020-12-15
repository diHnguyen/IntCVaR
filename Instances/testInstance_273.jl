edge = [1 7; 1 11; 1 13; 1 15; 2 3; 2 5; 2 6; 2 7; 2 9; 2 13; 2 15; 3 6; 3 10; 3 11; 3 13; 3 15; 4 8; 4 10; 4 11; 4 13; 4 15; 5 2; 5 3; 5 7; 5 10; 5 13; 5 14; 6 2; 6 3; 6 8; 6 10; 6 13; 7 3; 7 9; 7 10; 7 11; 7 12; 7 13; 8 3; 8 6; 8 11; 8 13; 8 14; 8 15; 9 3; 9 6; 9 7; 10 2; 10 7; 10 11; 11 5; 11 7; 11 13; 11 14; 12 2; 12 7; 12 8; 12 11; 13 3; 13 5; 13 8; 13 9; 13 14; 14 2; 14 5; 14 10; 14 11; 14 15]
cL_orig = [37.0, 450.0, 8.0, 4.0, 19.0, 121.0, 35.0, 170.0, 93.0, 167.0, 89.0, 94.0, 319.0, 52.0, 76.0, 167.0, 26.0, 61.0, 0.0, 65.0, 56.0, 49.0, 21.0, 66.0, 199.0, 137.0, 56.0, 76.0, 123.0, 67.0, 63.0, 113.0, 108.0, 52.0, 84.0, 177.0, 92.0, 15.0, 87.0, 21.0, 147.0, 101.0, 68.0, 84.0, 167.0, 23.0, 8.0, 88.0, 94.0, 27.0, 115.0, 127.0, 21.0, 83.0, 20.0, 179.0, 164.0, 0.0, 440.0, 328.0, 88.0, 8.0, 4.0, 537.0, 16.0, 127.0, 4.0, 40.0]
cU_orig = [89.0, 450.0, 8.0, 4.0, 19.0, 121.0, 35.0, 170.0, 509.0, 167.0, 349.0, 94.0, 319.0, 556.0, 76.0, 807.0, 124.0, 89.0, 602.0, 65.0, 56.0, 49.0, 177.0, 66.0, 199.0, 137.0, 670.0, 76.0, 123.0, 67.0, 283.0, 113.0, 108.0, 52.0, 84.0, 211.0, 92.0, 147.0, 177.0, 21.0, 147.0, 101.0, 192.0, 358.0, 167.0, 197.0, 144.0, 88.0, 94.0, 27.0, 115.0, 127.0, 21.0, 195.0, 20.0, 179.0, 164.0, 40.0, 440.0, 328.0, 88.0, 12.0, 4.0, 537.0, 280.0, 127.0, 30.0, 40.0]
d = [29.0, 37.0, 31.0, 43.0, 16.0, 37.0, 30.0, 5.0, 8.0, 36.0, 36.0, 7.0, 37.0, 35.0, 3.0, 46.0, 16.0, 28.0, 18.0, 44.0, 5.0, 46.0, 43.0, 8.0, 37.0, 10.0, 30.0, 25.0, 44.0, 9.0, 5.0, 47.0, 1.0, 43.0, 7.0, 40.0, 36.0, 31.0, 38.0, 25.0, 4.0, 5.0, 8.0, 4.0, 41.0, 47.0, 29.0, 18.0, 43.0, 15.0, 28.0, 30.0, 42.0, 20.0, 38.0, 6.0, 36.0, 32.0, 39.0, 27.0, 26.0, 4.0, 33.0, 35.0, 17.0, 38.0, 48.0, 33.0]
Len = length(d)

yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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