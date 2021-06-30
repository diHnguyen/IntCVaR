edge = [1 10; 1 12; 1 13; 1 14; 1 15; 1 16; 2 3; 2 8; 2 10; 2 11; 2 13; 2 14; 2 15; 3 7; 3 13; 3 15; 3 16; 3 17; 3 20; 4 5; 4 10; 5 2; 5 3; 5 4; 5 10; 5 13; 5 17; 5 19; 5 20; 6 4; 6 15; 6 16; 6 19; 7 3; 7 5; 7 15; 8 2; 8 4; 8 7; 8 10; 8 11; 8 12; 8 15; 8 17; 9 3; 9 6; 9 12; 9 13; 9 15; 9 17; 9 18; 9 19; 10 2; 10 4; 10 15; 10 20; 11 2; 11 4; 11 9; 11 15; 11 16; 11 20; 12 6; 12 11; 12 16; 12 17; 12 18; 13 4; 13 5; 13 9; 13 10; 13 15; 13 20; 14 4; 14 8; 14 10; 14 12; 14 15; 14 16; 14 19; 15 2; 15 3; 15 6; 15 9; 15 14; 15 19; 16 3; 16 5; 16 9; 16 17; 16 20; 17 2; 17 3; 17 8; 17 19; 18 4; 18 6; 18 7; 18 8; 18 11; 18 13; 18 17; 19 3; 19 7; 19 12; 19 14; 19 17]
cL_orig = [439.0, 540.0, 13.0, 483.0, 533.0, 571.0, 19.0, 4.0, 185.0, 384.0, 150.0, 19.0, 37.0, 11.0, 333.0, 582.0, 162.0, 592.0, 292.0, 7.0, 76.0, 27.0, 9.0, 33.0, 121.0, 244.0, 533.0, 268.0, 391.0, 14.0, 12.0, 106.0, 212.0, 8.0, 3.0, 38.0, 158.0, 169.0, 18.0, 12.0, 97.0, 183.0, 112.0, 345.0, 205.0, 45.0, 10.0, 32.0, 167.0, 318.0, 423.0, 346.0, 105.0, 181.0, 210.0, 372.0, 25.0, 222.0, 23.0, 84.0, 121.0, 35.0, 67.0, 25.0, 99.0, 85.0, 234.0, 126.0, 372.0, 136.0, 90.0, 2.0, 311.0, 210.0, 53.0, 117.0, 0.0, 0.0, 37.0, 9.0, 52.0, 171.0, 286.0, 245.0, 3.0, 121.0, 366.0, 97.0, 74.0, 29.0, 2.0, 562.0, 166.0, 382.0, 47.0, 355.0, 28.0, 390.0, 5.0, 34.0, 151.0, 50.0, 556.0, 225.0, 50.0, 244.0, 94.0]
cU_orig = [439.0, 540.0, 529.0, 483.0, 533.0, 571.0, 19.0, 4.0, 249.0, 384.0, 572.0, 627.0, 37.0, 203.0, 333.0, 582.0, 162.0, 592.0, 1334.0, 7.0, 76.0, 27.0, 23.0, 33.0, 121.0, 244.0, 533.0, 268.0, 391.0, 14.0, 24.0, 106.0, 212.0, 64.0, 3.0, 38.0, 158.0, 169.0, 18.0, 18.0, 97.0, 183.0, 112.0, 481.0, 259.0, 45.0, 34.0, 32.0, 167.0, 318.0, 423.0, 346.0, 105.0, 371.0, 210.0, 372.0, 633.0, 222.0, 23.0, 84.0, 121.0, 243.0, 67.0, 25.0, 99.0, 227.0, 234.0, 126.0, 372.0, 136.0, 90.0, 4.0, 311.0, 210.0, 53.0, 117.0, 2.0, 76.0, 37.0, 381.0, 52.0, 171.0, 286.0, 245.0, 41.0, 121.0, 366.0, 97.0, 74.0, 29.0, 372.0, 562.0, 848.0, 382.0, 47.0, 355.0, 842.0, 390.0, 5.0, 446.0, 151.0, 50.0, 556.0, 225.0, 50.0, 244.0, 94.0]
d = [47.0, 40.0, 19.0, 16.0, 35.0, 48.0, 50.0, 34.0, 10.0, 35.0, 22.0, 46.0, 41.0, 12.0, 18.0, 4.0, 28.0, 6.0, 14.0, 35.0, 24.0, 29.0, 19.0, 30.0, 35.0, 36.0, 18.0, 27.0, 29.0, 4.0, 47.0, 48.0, 35.0, 46.0, 19.0, 3.0, 16.0, 48.0, 43.0, 40.0, 19.0, 37.0, 15.0, 43.0, 31.0, 26.0, 21.0, 41.0, 8.0, 14.0, 38.0, 38.0, 33.0, 3.0, 27.0, 27.0, 16.0, 45.0, 22.0, 24.0, 7.0, 49.0, 23.0, 43.0, 48.0, 31.0, 50.0, 30.0, 17.0, 39.0, 45.0, 1.0, 20.0, 6.0, 12.0, 26.0, 50.0, 20.0, 22.0, 21.0, 40.0, 46.0, 9.0, 11.0, 1.0, 12.0, 7.0, 42.0, 23.0, 27.0, 18.0, 45.0, 16.0, 10.0, 28.0, 32.0, 17.0, 42.0, 5.0, 16.0, 41.0, 16.0, 2.0, 48.0, 19.0, 32.0, 36.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
