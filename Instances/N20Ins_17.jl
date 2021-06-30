edge = [1 10; 1 14; 1 15; 1 17; 1 19; 1 20; 2 3; 2 6; 2 10; 2 16; 2 18; 2 19; 3 6; 3 11; 3 14; 3 20; 4 2; 4 6; 4 7; 4 11; 4 13; 4 14; 4 19; 5 4; 5 7; 5 9; 5 11; 5 12; 5 14; 5 17; 5 19; 6 3; 6 4; 6 7; 6 8; 6 9; 6 10; 6 12; 6 13; 6 15; 6 16; 6 18; 7 4; 7 6; 7 8; 7 12; 7 14; 7 15; 7 16; 7 19; 8 5; 8 9; 8 11; 8 15; 9 7; 9 8; 9 18; 9 20; 10 2; 10 3; 10 9; 10 11; 10 15; 10 17; 10 20; 11 4; 11 8; 11 15; 11 16; 12 13; 12 14; 12 17; 13 4; 13 8; 13 9; 13 12; 13 15; 13 18; 13 20; 14 8; 14 10; 14 20; 15 4; 15 8; 15 17; 15 20; 16 2; 16 7; 16 18; 17 4; 17 5; 17 9; 17 15; 17 19; 17 20; 18 3; 18 5; 18 6; 18 7; 18 10; 18 11; 18 15; 18 16; 19 3; 19 4; 19 5; 19 11; 19 12; 19 15; 19 17; 19 18; 19 20]
cL_orig = [381.0, 318.0, 632.0, 700.0, 178.0, 622.0, 29.0, 63.0, 396.0, 128.0, 354.0, 701.0, 134.0, 120.0, 173.0, 696.0, 71.0, 9.0, 76.0, 254.0, 238.0, 395.0, 125.0, 20.0, 0.0, 148.0, 219.0, 4.0, 377.0, 63.0, 215.0, 60.0, 87.0, 27.0, 6.0, 98.0, 124.0, 46.0, 166.0, 213.0, 215.0, 539.0, 140.0, 34.0, 28.0, 213.0, 334.0, 80.0, 311.0, 71.0, 89.0, 41.0, 37.0, 15.0, 29.0, 19.0, 103.0, 32.0, 347.0, 251.0, 37.0, 19.0, 119.0, 211.0, 229.0, 157.0, 52.0, 105.0, 87.0, 50.0, 16.0, 173.0, 375.0, 45.0, 120.0, 26.0, 4.0, 152.0, 241.0, 300.0, 17.0, 20.0, 518.0, 12.0, 69.0, 126.0, 343.0, 284.0, 47.0, 148.0, 462.0, 150.0, 54.0, 63.0, 49.0, 746.0, 553.0, 271.0, 116.0, 248.0, 24.0, 77.0, 12.0, 418.0, 406.0, 519.0, 347.0, 98.0, 32.0, 1.0, 42.0, 37.0]
cU_orig = [381.0, 318.0, 632.0, 700.0, 348.0, 622.0, 29.0, 63.0, 396.0, 128.0, 354.0, 701.0, 134.0, 120.0, 173.0, 696.0, 71.0, 43.0, 122.0, 254.0, 298.0, 395.0, 125.0, 20.0, 128.0, 148.0, 219.0, 4.0, 377.0, 63.0, 753.0, 60.0, 87.0, 27.0, 14.0, 98.0, 124.0, 46.0, 166.0, 213.0, 215.0, 539.0, 140.0, 34.0, 28.0, 213.0, 334.0, 80.0, 311.0, 71.0, 89.0, 41.0, 37.0, 221.0, 29.0, 19.0, 103.0, 228.0, 347.0, 251.0, 37.0, 19.0, 119.0, 211.0, 229.0, 157.0, 52.0, 105.0, 87.0, 50.0, 16.0, 173.0, 375.0, 45.0, 120.0, 64.0, 60.0, 152.0, 241.0, 300.0, 17.0, 20.0, 518.0, 172.0, 69.0, 126.0, 343.0, 284.0, 47.0, 372.0, 462.0, 150.0, 54.0, 63.0, 49.0, 746.0, 553.0, 271.0, 146.0, 284.0, 24.0, 77.0, 12.0, 418.0, 406.0, 519.0, 347.0, 98.0, 32.0, 1.0, 42.0, 37.0]
d = [3.0, 39.0, 4.0, 9.0, 23.0, 35.0, 10.0, 18.0, 44.0, 25.0, 40.0, 30.0, 21.0, 13.0, 17.0, 35.0, 46.0, 6.0, 10.0, 33.0, 42.0, 12.0, 40.0, 8.0, 40.0, 13.0, 2.0, 45.0, 24.0, 23.0, 6.0, 21.0, 33.0, 16.0, 14.0, 34.0, 27.0, 28.0, 20.0, 11.0, 39.0, 20.0, 32.0, 16.0, 19.0, 14.0, 46.0, 44.0, 6.0, 44.0, 36.0, 43.0, 28.0, 13.0, 40.0, 9.0, 22.0, 17.0, 12.0, 11.0, 3.0, 30.0, 10.0, 43.0, 2.0, 37.0, 1.0, 6.0, 43.0, 35.0, 8.0, 43.0, 41.0, 25.0, 38.0, 21.0, 8.0, 19.0, 24.0, 25.0, 48.0, 19.0, 15.0, 25.0, 46.0, 31.0, 28.0, 24.0, 31.0, 17.0, 12.0, 50.0, 19.0, 21.0, 20.0, 46.0, 35.0, 2.0, 8.0, 50.0, 8.0, 24.0, 45.0, 40.0, 34.0, 28.0, 48.0, 21.0, 20.0, 41.0, 9.0, 1.0]
Len = length(d)

yy = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

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
