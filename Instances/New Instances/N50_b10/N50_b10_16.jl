edge = [1 3; 1 18; 1 27; 1 28; 1 47; 2 4; 2 5; 2 16; 2 31; 2 45; 3 8; 3 25; 3 33; 3 37; 3 39; 4 5; 4 14; 4 20; 4 43; 5 11; 5 28; 5 41; 5 50; 6 16; 6 21; 6 40; 7 5; 7 9; 7 16; 7 20; 7 23; 7 32; 7 34; 7 43; 8 6; 8 11; 8 16; 8 18; 8 29; 9 7; 9 21; 9 22; 9 37; 9 43; 10 20; 10 23; 10 43; 10 45; 11 3; 11 7; 11 15; 11 37; 11 42; 12 4; 12 13; 12 19; 12 20; 12 22; 12 29; 12 31; 12 40; 13 21; 13 31; 14 6; 14 13; 14 21; 14 40; 14 47; 15 8; 15 13; 15 18; 15 27; 15 31; 15 41; 15 46; 15 47; 16 4; 16 6; 16 11; 16 37; 17 3; 17 5; 17 7; 17 22; 17 24; 17 35; 17 41; 17 48; 18 49; 19 6; 19 16; 19 24; 19 48; 20 30; 20 38; 20 40; 20 44; 20 46; 20 47; 21 3; 21 13; 21 27; 22 19; 22 32; 23 35; 23 40; 24 27; 24 34; 24 45; 25 2; 25 3; 25 8; 25 12; 25 40; 25 50; 26 5; 26 10; 26 14; 26 34; 26 35; 26 38; 27 7; 27 8; 27 11; 27 16; 27 22; 27 42; 27 44; 28 4; 28 11; 28 12; 28 18; 28 25; 28 33; 28 48; 29 5; 29 11; 29 17; 29 28; 29 38; 29 39; 29 46; 30 17; 30 26; 30 47; 30 48; 31 4; 31 9; 31 15; 31 24; 32 5; 32 7; 32 29; 32 43; 33 3; 33 26; 33 30; 33 31; 33 40; 33 41; 33 44; 34 12; 34 19; 34 25; 34 26; 34 31; 35 15; 35 17; 36 6; 36 9; 36 10; 36 12; 36 18; 36 26; 36 37; 36 38; 36 39; 37 9; 37 12; 37 27; 37 28; 37 35; 37 36; 37 38; 37 44; 37 50; 38 11; 38 20; 38 28; 39 5; 39 6; 39 10; 39 12; 39 16; 39 31; 39 35; 39 42; 40 5; 40 30; 40 46; 41 7; 41 16; 41 23; 41 24; 41 29; 41 31; 42 21; 42 38; 42 39; 42 40; 42 50; 43 4; 43 22; 43 49; 44 6; 44 14; 44 17; 44 20; 45 5; 45 25; 45 39; 45 41; 45 49; 45 50; 46 5; 46 6; 46 10; 46 19; 46 24; 46 37; 46 40; 46 42; 46 44; 46 45; 47 3; 47 16; 47 17; 47 31; 48 7; 48 15; 48 20; 48 24; 48 26; 48 27; 48 33; 49 3; 49 9]
cL_orig = [20.0, 156.0, 255.0, 270.0, 456.0, 25.0, 30.0, 137.0, 285.0, 433.0, 42.0, 206.0, 304.0, 332.0, 362.0, 10.0, 104.0, 158.0, 383.0, 64.0, 220.0, 364.0, 442.0, 101.0, 146.0, 344.0, 15.0, 24.0, 94.0, 119.0, 163.0, 242.0, 274.0, 362.0, 15.0, 31.0, 83.0, 99.0, 211.0, 15.0, 120.0, 120.0, 277.0, 345.0, 100.0, 115.0, 325.0, 345.0, 78.0, 40.0, 35.0, 255.0, 308.0, 77.0, 2.0, 68.0, 67.0, 101.0, 170.0, 185.0, 282.0, 84.0, 177.0, 82.0, 8.0, 67.0, 258.0, 331.0, 68.0, 20.0, 29.0, 122.0, 165.0, 248.0, 311.0, 319.0, 125.0, 103.0, 50.0, 211.0, 137.0, 124.0, 102.0, 51.0, 67.0, 183.0, 237.0, 313.0, 308.0, 127.0, 19.0, 49.0, 295.0, 105.0, 175.0, 203.0, 245.0, 261.0, 267.0, 180.0, 71.0, 65.0, 26.0, 104.0, 114.0, 175.0, 35.0, 103.0, 211.0, 227.0, 224.0, 159.0, 130.0, 137.0, 252.0, 209.0, 164.0, 125.0, 72.0, 94.0, 121.0, 202.0, 191.0, 160.0, 107.0, 41.0, 146.0, 165.0, 235.0, 168.0, 159.0, 97.0, 35.0, 47.0, 203.0, 235.0, 180.0, 121.0, 2.0, 92.0, 91.0, 168.0, 127.0, 41.0, 168.0, 177.0, 274.0, 216.0, 164.0, 57.0, 273.0, 239.0, 25.0, 111.0, 302.0, 66.0, 34.0, 17.0, 72.0, 84.0, 107.0, 225.0, 152.0, 94.0, 82.0, 31.0, 203.0, 185.0, 304.0, 265.0, 256.0, 235.0, 171.0, 97.0, 14.0, 19.0, 25.0, 278.0, 250.0, 103.0, 89.0, 19.0, 13.0, 9.0, 55.0, 132.0, 265.0, 183.0, 102.0, 336.0, 334.0, 287.0, 269.0, 227.0, 77.0, 44.0, 30.0, 355.0, 103.0, 61.0, 337.0, 246.0, 179.0, 172.0, 119.0, 95.0, 207.0, 37.0, 27.0, 23.0, 66.0, 385.0, 212.0, 56.0, 376.0, 300.0, 263.0, 238.0, 405.0, 195.0, 62.0, 40.0, 42.0, 49.0, 405.0, 395.0, 360.0, 270.0, 216.0, 90.0, 64.0, 43.0, 19.0, 0.0, 427.0, 309.0, 304.0, 151.0, 399.0, 321.0, 275.0, 237.0, 209.0, 202.0, 148.0, 455.0, 400.0]
cU_orig = [20.0, 176.0, 255.0, 270.0, 456.0, 25.0, 30.0, 137.0, 297.0, 433.0, 58.0, 226.0, 304.0, 346.0, 362.0, 10.0, 104.0, 158.0, 393.0, 64.0, 240.0, 364.0, 460.0, 101.0, 146.0, 344.0, 15.0, 24.0, 94.0, 131.0, 163.0, 254.0, 274.0, 362.0, 15.0, 31.0, 83.0, 99.0, 211.0, 29.0, 120.0, 132.0, 277.0, 345.0, 100.0, 135.0, 325.0, 345.0, 78.0, 40.0, 35.0, 255.0, 308.0, 77.0, 14.0, 68.0, 85.0, 101.0, 170.0, 195.0, 282.0, 84.0, 177.0, 82.0, 8.0, 67.0, 258.0, 331.0, 68.0, 20.0, 29.0, 122.0, 165.0, 268.0, 311.0, 319.0, 125.0, 103.0, 50.0, 211.0, 137.0, 124.0, 102.0, 51.0, 67.0, 183.0, 237.0, 313.0, 308.0, 127.0, 33.0, 49.0, 295.0, 105.0, 193.0, 203.0, 245.0, 261.0, 267.0, 180.0, 87.0, 65.0, 26.0, 104.0, 134.0, 175.0, 35.0, 103.0, 211.0, 227.0, 224.0, 179.0, 130.0, 157.0, 252.0, 209.0, 164.0, 125.0, 90.0, 94.0, 121.0, 202.0, 191.0, 160.0, 107.0, 55.0, 146.0, 165.0, 235.0, 168.0, 159.0, 97.0, 35.0, 47.0, 203.0, 235.0, 180.0, 121.0, 16.0, 92.0, 105.0, 168.0, 127.0, 41.0, 168.0, 193.0, 274.0, 216.0, 164.0, 73.0, 273.0, 253.0, 45.0, 111.0, 302.0, 66.0, 34.0, 27.0, 72.0, 84.0, 107.0, 225.0, 152.0, 94.0, 82.0, 31.0, 203.0, 185.0, 304.0, 265.0, 256.0, 235.0, 185.0, 97.0, 14.0, 19.0, 25.0, 278.0, 250.0, 103.0, 89.0, 19.0, 13.0, 19.0, 75.0, 132.0, 265.0, 183.0, 102.0, 336.0, 334.0, 287.0, 269.0, 227.0, 77.0, 44.0, 30.0, 355.0, 103.0, 61.0, 337.0, 246.0, 179.0, 172.0, 119.0, 95.0, 207.0, 51.0, 27.0, 23.0, 84.0, 385.0, 212.0, 56.0, 376.0, 300.0, 275.0, 238.0, 405.0, 211.0, 62.0, 40.0, 42.0, 49.0, 405.0, 415.0, 360.0, 270.0, 216.0, 90.0, 64.0, 43.0, 31.0, 18.0, 447.0, 309.0, 304.0, 167.0, 415.0, 341.0, 275.0, 237.0, 223.0, 218.0, 148.0, 469.0, 400.0]
d = [9.0, 20.0, 7.0, 20.0, 13.0, 16.0, 12.0, 16.0, 16.0, 7.0, 1.0, 2.0, 16.0, 12.0, 6.0, 18.0, 3.0, 9.0, 8.0, 15.0, 20.0, 6.0, 19.0, 11.0, 16.0, 6.0, 19.0, 11.0, 5.0, 17.0, 17.0, 12.0, 12.0, 10.0, 6.0, 3.0, 18.0, 2.0, 10.0, 17.0, 10.0, 20.0, 18.0, 5.0, 4.0, 19.0, 20.0, 11.0, 18.0, 14.0, 11.0, 12.0, 5.0, 18.0, 18.0, 15.0, 19.0, 9.0, 15.0, 9.0, 1.0, 19.0, 16.0, 16.0, 11.0, 8.0, 17.0, 9.0, 10.0, 4.0, 8.0, 4.0, 9.0, 20.0, 6.0, 18.0, 16.0, 14.0, 14.0, 17.0, 2.0, 5.0, 7.0, 5.0, 17.0, 3.0, 13.0, 19.0, 2.0, 4.0, 15.0, 14.0, 11.0, 1.0, 18.0, 13.0, 19.0, 7.0, 7.0, 18.0, 13.0, 4.0, 8.0, 11.0, 9.0, 13.0, 13.0, 13.0, 12.0, 2.0, 6.0, 14.0, 11.0, 2.0, 19.0, 2.0, 10.0, 11.0, 19.0, 16.0, 13.0, 2.0, 14.0, 1.0, 6.0, 10.0, 16.0, 16.0, 11.0, 8.0, 6.0, 12.0, 15.0, 18.0, 13.0, 17.0, 19.0, 19.0, 3.0, 6.0, 13.0, 19.0, 7.0, 20.0, 16.0, 3.0, 6.0, 7.0, 18.0, 11.0, 15.0, 18.0, 6.0, 1.0, 19.0, 4.0, 18.0, 14.0, 17.0, 17.0, 19.0, 4.0, 6.0, 14.0, 19.0, 17.0, 1.0, 10.0, 5.0, 15.0, 4.0, 12.0, 5.0, 13.0, 11.0, 15.0, 14.0, 13.0, 14.0, 8.0, 18.0, 5.0, 11.0, 4.0, 6.0, 20.0, 1.0, 19.0, 3.0, 1.0, 4.0, 6.0, 4.0, 12.0, 17.0, 18.0, 16.0, 5.0, 19.0, 16.0, 5.0, 16.0, 3.0, 14.0, 10.0, 10.0, 20.0, 6.0, 9.0, 11.0, 4.0, 2.0, 16.0, 5.0, 14.0, 8.0, 14.0, 19.0, 18.0, 5.0, 19.0, 9.0, 1.0, 18.0, 18.0, 14.0, 5.0, 18.0, 9.0, 20.0, 2.0, 15.0, 7.0, 11.0, 13.0, 17.0, 20.0, 7.0, 18.0, 12.0, 17.0, 13.0, 17.0, 17.0, 15.0, 5.0, 7.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
b = 10
last_node = maximum(edge)
β = 0.323