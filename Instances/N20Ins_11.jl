edge = [1 10; 1 12; 1 14; 1 15; 1 16; 1 19; 2 9; 2 10; 2 13; 2 20; 3 5; 3 7; 3 8; 3 10; 3 12; 3 13; 3 17; 3 18; 3 20; 4 2; 4 3; 4 5; 4 12; 4 18; 5 4; 5 8; 5 10; 5 14; 5 15; 5 18; 5 20; 6 5; 6 7; 6 16; 6 17; 6 18; 6 19; 7 2; 7 6; 7 12; 7 14; 8 3; 8 5; 8 12; 8 13; 8 15; 8 17; 8 20; 9 2; 9 8; 9 10; 9 14; 9 19; 10 2; 10 11; 10 13; 10 19; 11 4; 11 5; 11 6; 11 8; 11 12; 11 13; 11 15; 11 19; 12 3; 12 4; 12 13; 12 15; 12 16; 12 17; 12 19; 13 3; 13 4; 13 18; 14 5; 14 6; 14 7; 14 9; 14 12; 14 15; 14 17; 14 18; 15 5; 15 7; 15 10; 15 11; 15 13; 15 20; 16 8; 16 9; 16 10; 16 12; 16 19; 16 20; 17 4; 17 13; 17 14; 17 20; 18 3; 18 9; 18 16; 18 19; 18 20; 19 2; 19 6; 19 8; 19 9; 19 15; 19 16; 19 17; 19 20]
cL_orig = [10.0, 413.0, 16.0, 189.0, 479.0, 576.0, 338.0, 5.0, 138.0, 361.0, 36.0, 100.0, 30.0, 48.0, 324.0, 411.0, 51.0, 0.0, 31.0, 30.0, 12.0, 19.0, 130.0, 119.0, 9.0, 145.0, 13.0, 333.0, 2.0, 391.0, 125.0, 5.0, 20.0, 428.0, 379.0, 311.0, 147.0, 97.0, 22.0, 115.0, 286.0, 232.0, 131.0, 10.0, 238.0, 190.0, 358.0, 146.0, 216.0, 34.0, 50.0, 53.0, 308.0, 319.0, 20.0, 14.0, 19.0, 82.0, 138.0, 178.0, 80.0, 39.0, 69.0, 126.0, 335.0, 428.0, 38.0, 0.0, 51.0, 119.0, 214.0, 320.0, 40.0, 428.0, 51.0, 117.0, 374.0, 208.0, 231.0, 57.0, 5.0, 103.0, 122.0, 53.0, 196.0, 91.0, 3.0, 31.0, 52.0, 270.0, 242.0, 57.0, 1.0, 88.0, 95.0, 550.0, 110.0, 35.0, 119.0, 338.0, 124.0, 64.0, 36.0, 92.0, 720.0, 326.0, 134.0, 292.0, 87.0, 6.0, 12.0, 5.0]
cU_orig = [10.0, 413.0, 16.0, 189.0, 479.0, 576.0, 338.0, 485.0, 138.0, 787.0, 36.0, 100.0, 48.0, 316.0, 324.0, 411.0, 331.0, 34.0, 57.0, 136.0, 12.0, 63.0, 130.0, 119.0, 37.0, 145.0, 39.0, 333.0, 18.0, 391.0, 125.0, 43.0, 20.0, 428.0, 379.0, 311.0, 147.0, 97.0, 22.0, 115.0, 286.0, 232.0, 131.0, 162.0, 238.0, 190.0, 358.0, 146.0, 216.0, 34.0, 50.0, 139.0, 308.0, 319.0, 20.0, 14.0, 19.0, 270.0, 138.0, 178.0, 80.0, 39.0, 69.0, 126.0, 335.0, 428.0, 38.0, 52.0, 125.0, 119.0, 214.0, 320.0, 40.0, 428.0, 51.0, 147.0, 374.0, 208.0, 231.0, 57.0, 53.0, 103.0, 122.0, 53.0, 196.0, 143.0, 83.0, 151.0, 52.0, 270.0, 242.0, 313.0, 5.0, 88.0, 95.0, 550.0, 182.0, 35.0, 119.0, 338.0, 124.0, 64.0, 36.0, 92.0, 720.0, 326.0, 148.0, 292.0, 87.0, 268.0, 70.0, 5.0]
d = [42.0, 30.0, 21.0, 17.0, 28.0, 5.0, 28.0, 47.0, 50.0, 18.0, 23.0, 12.0, 27.0, 39.0, 27.0, 45.0, 7.0, 31.0, 12.0, 11.0, 34.0, 49.0, 47.0, 32.0, 16.0, 12.0, 42.0, 45.0, 33.0, 37.0, 6.0, 30.0, 1.0, 37.0, 14.0, 9.0, 45.0, 33.0, 15.0, 16.0, 45.0, 15.0, 48.0, 30.0, 12.0, 45.0, 24.0, 4.0, 22.0, 19.0, 13.0, 15.0, 11.0, 35.0, 12.0, 45.0, 36.0, 29.0, 39.0, 16.0, 36.0, 48.0, 40.0, 16.0, 35.0, 17.0, 29.0, 3.0, 8.0, 37.0, 35.0, 45.0, 16.0, 44.0, 34.0, 42.0, 47.0, 33.0, 20.0, 18.0, 12.0, 49.0, 45.0, 21.0, 14.0, 11.0, 14.0, 43.0, 49.0, 37.0, 33.0, 44.0, 41.0, 13.0, 31.0, 2.0, 36.0, 33.0, 6.0, 47.0, 14.0, 38.0, 46.0, 49.0, 8.0, 21.0, 21.0, 14.0, 8.0, 48.0, 37.0, 44.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
delta2 = 2
last_node = maximum(edge)
