edge = [1 4; 1 5; 1 7; 1 8; 1 12; 1 19; 1 20; 2 4; 2 5; 2 12; 3 6; 3 9; 3 11; 3 19; 4 2; 4 7; 4 8; 4 12; 4 15; 4 19; 5 2; 5 3; 5 7; 5 14; 5 15; 5 19; 6 3; 6 4; 6 19; 6 20; 7 3; 7 10; 7 12; 7 15; 7 16; 8 6; 8 10; 8 13; 9 4; 9 10; 9 16; 9 17; 9 18; 9 20; 10 2; 10 5; 10 6; 10 7; 10 9; 10 13; 10 16; 10 17; 10 19; 11 2; 11 4; 11 14; 11 15; 11 17; 11 18; 12 4; 12 7; 12 10; 12 11; 12 13; 12 15; 13 2; 13 8; 13 14; 13 15; 13 17; 14 3; 14 6; 14 9; 14 12; 14 13; 14 17; 14 18; 15 3; 15 4; 15 8; 15 10; 15 11; 15 12; 15 14; 15 19; 16 4; 16 5; 16 12; 16 13; 16 14; 16 18; 16 19; 17 11; 17 12; 17 18; 17 20; 18 2; 18 3; 18 4; 18 16; 18 17; 18 19; 19 2; 19 3; 19 7; 19 10; 19 12; 19 15]
cL_orig = [37.0, 163.0, 76.0, 252.0, 329.0, 106.0, 135.0, 8.0, 10.0, 0.0, 85.0, 212.0, 135.0, 134.0, 9.0, 39.0, 37.0, 153.0, 73.0, 160.0, 123.0, 0.0, 49.0, 291.0, 138.0, 540.0, 119.0, 9.0, 378.0, 259.0, 184.0, 113.0, 202.0, 46.0, 411.0, 16.0, 5.0, 33.0, 157.0, 31.0, 245.0, 244.0, 160.0, 89.0, 51.0, 40.0, 98.0, 105.0, 11.0, 107.0, 77.0, 57.0, 217.0, 376.0, 283.0, 67.0, 195.0, 202.0, 263.0, 69.0, 163.0, 67.0, 21.0, 25.0, 31.0, 472.0, 85.0, 3.0, 8.0, 72.0, 106.0, 290.0, 204.0, 0.0, 8.0, 105.0, 73.0, 158.0, 528.0, 288.0, 122.0, 62.0, 31.0, 45.0, 74.0, 491.0, 404.0, 96.0, 5.0, 89.0, 89.0, 96.0, 284.0, 2.0, 31.0, 9.0, 640.0, 700.0, 658.0, 10.0, 23.0, 20.0, 258.0, 150.0, 196.0, 261.0, 191.0, 30.0]
cU_orig = [37.0, 163.0, 80.0, 252.0, 329.0, 106.0, 135.0, 8.0, 120.0, 6.0, 95.0, 212.0, 135.0, 134.0, 175.0, 39.0, 37.0, 153.0, 73.0, 160.0, 123.0, 22.0, 49.0, 291.0, 138.0, 540.0, 119.0, 9.0, 378.0, 259.0, 184.0, 113.0, 202.0, 46.0, 411.0, 16.0, 5.0, 33.0, 157.0, 31.0, 245.0, 244.0, 160.0, 89.0, 51.0, 276.0, 98.0, 105.0, 11.0, 107.0, 77.0, 57.0, 217.0, 376.0, 283.0, 67.0, 195.0, 202.0, 263.0, 131.0, 163.0, 67.0, 21.0, 25.0, 31.0, 472.0, 85.0, 9.0, 94.0, 72.0, 670.0, 290.0, 204.0, 178.0, 14.0, 105.0, 73.0, 158.0, 528.0, 288.0, 122.0, 62.0, 31.0, 45.0, 74.0, 491.0, 404.0, 96.0, 53.0, 89.0, 89.0, 96.0, 284.0, 92.0, 31.0, 117.0, 640.0, 700.0, 658.0, 10.0, 23.0, 20.0, 368.0, 150.0, 196.0, 261.0, 191.0, 74.0]
d = [13.0, 8.0, 47.0, 22.0, 19.0, 24.0, 2.0, 17.0, 8.0, 13.0, 27.0, 32.0, 32.0, 42.0, 3.0, 4.0, 18.0, 42.0, 16.0, 45.0, 7.0, 1.0, 8.0, 37.0, 50.0, 42.0, 31.0, 25.0, 43.0, 1.0, 19.0, 34.0, 28.0, 19.0, 34.0, 12.0, 18.0, 14.0, 44.0, 21.0, 20.0, 13.0, 6.0, 21.0, 8.0, 17.0, 2.0, 21.0, 25.0, 24.0, 23.0, 36.0, 19.0, 16.0, 14.0, 35.0, 19.0, 35.0, 6.0, 24.0, 10.0, 44.0, 47.0, 28.0, 29.0, 48.0, 19.0, 4.0, 30.0, 10.0, 39.0, 11.0, 5.0, 1.0, 35.0, 38.0, 49.0, 21.0, 18.0, 9.0, 44.0, 11.0, 46.0, 7.0, 46.0, 6.0, 6.0, 25.0, 10.0, 41.0, 36.0, 41.0, 6.0, 26.0, 30.0, 28.0, 24.0, 21.0, 2.0, 3.0, 25.0, 26.0, 33.0, 35.0, 5.0, 26.0, 1.0, 27.0]
Len = length(d)

yy = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =20

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
