edge = [1 4; 1 6; 1 41; 1 47; 2 9; 2 39; 3 10; 3 13; 3 29; 3 47; 4 7; 4 8; 4 27; 5 2; 5 3; 5 11; 5 18; 5 32; 6 3; 6 9; 6 17; 6 19; 6 30; 6 34; 6 36; 6 40; 7 13; 7 23; 8 15; 8 19; 8 21; 8 30; 8 34; 8 38; 8 48; 8 49; 9 12; 9 20; 9 27; 9 34; 9 47; 10 12; 11 6; 11 8; 11 15; 12 2; 12 24; 12 37; 12 47; 13 2; 13 3; 13 6; 13 11; 13 33; 13 36; 13 45; 14 3; 14 8; 14 26; 14 36; 14 37; 14 38; 14 43; 14 44; 15 11; 15 21; 15 32; 15 42; 16 3; 16 4; 16 5; 16 21; 16 29; 17 20; 17 27; 17 34; 17 49; 17 50; 18 8; 18 38; 18 50; 19 29; 19 37; 19 49; 20 30; 20 40; 20 46; 20 47; 21 8; 21 20; 21 23; 21 37; 22 6; 22 11; 22 46; 23 9; 23 32; 23 37; 23 44; 24 23; 24 44; 24 47; 25 6; 25 29; 25 31; 25 35; 25 36; 25 40; 26 17; 26 30; 26 33; 26 43; 27 2; 27 3; 27 9; 27 13; 27 16; 27 38; 28 25; 28 34; 28 37; 28 47; 28 50; 29 6; 29 14; 29 15; 29 23; 29 25; 29 34; 29 38; 29 40; 29 41; 30 8; 30 19; 30 24; 30 31; 30 36; 31 36; 31 44; 31 48; 32 3; 32 9; 33 9; 33 12; 33 23; 33 24; 33 30; 33 35; 34 13; 34 20; 34 33; 34 35; 35 2; 35 14; 35 29; 35 38; 36 7; 36 8; 36 35; 36 38; 36 39; 37 9; 37 14; 37 24; 37 25; 37 29; 37 38; 38 41; 38 48; 39 4; 39 13; 39 14; 39 19; 39 25; 39 41; 39 46; 40 17; 40 21; 40 24; 40 32; 40 35; 41 8; 41 22; 42 12; 42 15; 42 38; 42 47; 43 2; 43 5; 43 7; 43 22; 43 29; 43 38; 43 40; 43 41; 43 47; 44 21; 44 29; 44 30; 44 40; 45 3; 45 6; 45 12; 45 32; 45 34; 45 35; 46 14; 46 28; 47 6; 47 12; 47 17; 47 22; 47 37; 47 40; 47 45; 47 50; 48 4; 48 39; 48 41; 48 43; 48 49; 49 10; 49 11; 49 12; 49 35; 49 36; 49 40; 49 42]
cL_orig = [30.0, 39.0, 399.0, 458.0, 72.0, 363.0, 61.0, 91.0, 251.0, 431.0, 29.0, 39.0, 232.0, 33.0, 7.0, 56.0, 125.0, 275.0, 32.0, 27.0, 110.0, 121.0, 231.0, 278.0, 305.0, 331.0, 50.0, 163.0, 65.0, 110.0, 132.0, 224.0, 252.0, 299.0, 403.0, 407.0, 34.0, 107.0, 182.0, 248.0, 377.0, 19.0, 50.0, 28.0, 41.0, 102.0, 111.0, 245.0, 352.0, 109.0, 91.0, 60.0, 24.0, 202.0, 229.0, 309.0, 103.0, 56.0, 125.0, 225.0, 228.0, 235.0, 291.0, 298.0, 25.0, 53.0, 160.0, 268.0, 127.0, 120.0, 96.0, 47.0, 120.0, 19.0, 99.0, 173.0, 318.0, 328.0, 96.0, 201.0, 321.0, 92.0, 184.0, 292.0, 90.0, 187.0, 255.0, 264.0, 128.0, 7.0, 11.0, 152.0, 157.0, 107.0, 239.0, 141.0, 86.0, 144.0, 207.0, 7.0, 201.0, 221.0, 176.0, 31.0, 54.0, 102.0, 102.0, 142.0, 93.0, 44.0, 74.0, 171.0, 244.0, 238.0, 185.0, 130.0, 107.0, 103.0, 25.0, 56.0, 84.0, 188.0, 218.0, 231.0, 145.0, 140.0, 59.0, 35.0, 55.0, 80.0, 107.0, 117.0, 209.0, 109.0, 53.0, 8.0, 63.0, 42.0, 127.0, 167.0, 293.0, 235.0, 237.0, 211.0, 95.0, 88.0, 25.0, 16.0, 208.0, 140.0, 6.0, 4.0, 335.0, 215.0, 52.0, 30.0, 292.0, 275.0, 1.0, 19.0, 33.0, 275.0, 235.0, 126.0, 124.0, 75.0, 5.0, 22.0, 105.0, 349.0, 254.0, 252.0, 201.0, 140.0, 16.0, 73.0, 225.0, 182.0, 153.0, 78.0, 41.0, 319.0, 177.0, 302.0, 258.0, 44.0, 55.0, 399.0, 375.0, 350.0, 214.0, 137.0, 46.0, 17.0, 24.0, 39.0, 226.0, 139.0, 137.0, 39.0, 422.0, 375.0, 323.0, 119.0, 109.0, 88.0, 319.0, 185.0, 403.0, 350.0, 299.0, 243.0, 93.0, 74.0, 11.0, 20.0, 437.0, 88.0, 64.0, 42.0, 8.0, 395.0, 385.0, 375.0, 127.0, 133.0, 86.0, 73.0]
cU_orig = [30.0, 53.0, 399.0, 470.0, 72.0, 379.0, 71.0, 111.0, 261.0, 443.0, 39.0, 39.0, 232.0, 33.0, 27.0, 74.0, 125.0, 275.0, 32.0, 41.0, 110.0, 131.0, 247.0, 278.0, 305.0, 349.0, 66.0, 163.0, 65.0, 110.0, 132.0, 224.0, 270.0, 299.0, 403.0, 407.0, 34.0, 107.0, 182.0, 248.0, 377.0, 19.0, 50.0, 28.0, 41.0, 102.0, 131.0, 245.0, 352.0, 119.0, 103.0, 72.0, 24.0, 202.0, 229.0, 327.0, 115.0, 56.0, 125.0, 225.0, 228.0, 253.0, 291.0, 298.0, 45.0, 67.0, 172.0, 280.0, 127.0, 120.0, 116.0, 61.0, 140.0, 37.0, 99.0, 173.0, 318.0, 328.0, 96.0, 201.0, 321.0, 104.0, 184.0, 312.0, 102.0, 205.0, 255.0, 276.0, 128.0, 7.0, 27.0, 170.0, 157.0, 107.0, 239.0, 141.0, 86.0, 144.0, 221.0, 7.0, 201.0, 231.0, 194.0, 41.0, 70.0, 102.0, 112.0, 160.0, 93.0, 44.0, 74.0, 171.0, 264.0, 238.0, 185.0, 150.0, 107.0, 115.0, 43.0, 56.0, 98.0, 188.0, 218.0, 231.0, 145.0, 140.0, 59.0, 45.0, 55.0, 94.0, 117.0, 117.0, 225.0, 109.0, 67.0, 8.0, 63.0, 60.0, 127.0, 177.0, 293.0, 235.0, 237.0, 211.0, 95.0, 88.0, 37.0, 16.0, 208.0, 140.0, 24.0, 16.0, 335.0, 215.0, 70.0, 30.0, 292.0, 275.0, 19.0, 19.0, 33.0, 289.0, 235.0, 142.0, 124.0, 75.0, 25.0, 36.0, 105.0, 349.0, 266.0, 252.0, 201.0, 140.0, 16.0, 73.0, 225.0, 202.0, 163.0, 78.0, 53.0, 335.0, 195.0, 302.0, 276.0, 44.0, 55.0, 413.0, 375.0, 362.0, 214.0, 137.0, 46.0, 37.0, 24.0, 49.0, 236.0, 151.0, 137.0, 39.0, 422.0, 395.0, 339.0, 139.0, 109.0, 106.0, 319.0, 185.0, 413.0, 350.0, 299.0, 263.0, 103.0, 74.0, 27.0, 32.0, 437.0, 88.0, 80.0, 62.0, 8.0, 395.0, 385.0, 375.0, 145.0, 133.0, 86.0, 73.0]
d = [3.0, 9.0, 20.0, 13.0, 4.0, 16.0, 5.0, 18.0, 12.0, 15.0, 19.0, 15.0, 6.0, 7.0, 14.0, 10.0, 12.0, 15.0, 15.0, 8.0, 9.0, 11.0, 19.0, 4.0, 12.0, 10.0, 6.0, 13.0, 17.0, 20.0, 7.0, 9.0, 7.0, 4.0, 20.0, 6.0, 8.0, 11.0, 10.0, 9.0, 16.0, 14.0, 11.0, 20.0, 7.0, 9.0, 11.0, 13.0, 5.0, 6.0, 1.0, 18.0, 14.0, 13.0, 16.0, 13.0, 10.0, 4.0, 7.0, 12.0, 14.0, 6.0, 18.0, 8.0, 9.0, 10.0, 4.0, 18.0, 7.0, 7.0, 2.0, 9.0, 14.0, 18.0, 15.0, 4.0, 20.0, 6.0, 5.0, 14.0, 10.0, 16.0, 7.0, 11.0, 12.0, 6.0, 4.0, 19.0, 14.0, 4.0, 14.0, 17.0, 19.0, 15.0, 17.0, 17.0, 19.0, 20.0, 16.0, 3.0, 18.0, 2.0, 11.0, 11.0, 1.0, 6.0, 5.0, 2.0, 10.0, 6.0, 15.0, 17.0, 5.0, 17.0, 6.0, 20.0, 9.0, 9.0, 20.0, 5.0, 6.0, 20.0, 7.0, 11.0, 6.0, 8.0, 1.0, 2.0, 1.0, 13.0, 2.0, 3.0, 20.0, 8.0, 11.0, 7.0, 3.0, 3.0, 12.0, 9.0, 19.0, 12.0, 5.0, 13.0, 4.0, 15.0, 11.0, 10.0, 11.0, 8.0, 16.0, 12.0, 14.0, 19.0, 1.0, 1.0, 9.0, 19.0, 7.0, 14.0, 1.0, 14.0, 20.0, 20.0, 3.0, 19.0, 9.0, 19.0, 12.0, 1.0, 20.0, 16.0, 3.0, 19.0, 10.0, 17.0, 9.0, 2.0, 20.0, 4.0, 11.0, 2.0, 4.0, 15.0, 11.0, 4.0, 7.0, 20.0, 14.0, 16.0, 5.0, 1.0, 20.0, 5.0, 5.0, 5.0, 8.0, 16.0, 6.0, 19.0, 3.0, 19.0, 7.0, 20.0, 16.0, 10.0, 3.0, 2.0, 19.0, 1.0, 1.0, 20.0, 19.0, 17.0, 2.0, 3.0, 17.0, 13.0, 8.0, 4.0, 20.0, 11.0, 11.0, 8.0, 5.0, 14.0, 12.0, 4.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 50
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 2.0
b = 7
last_node = maximum(edge)
β = 0.604
