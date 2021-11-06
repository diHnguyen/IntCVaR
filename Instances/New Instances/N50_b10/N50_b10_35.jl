edge = [1 3; 1 7; 1 8; 1 16; 1 26; 1 31; 1 33; 2 6; 2 9; 2 11; 2 18; 2 36; 2 46; 3 19; 3 22; 3 27; 3 29; 3 39; 3 44; 4 27; 4 32; 4 39; 4 40; 5 8; 5 20; 5 41; 5 42; 5 49; 5 50; 6 18; 6 45; 6 46; 7 18; 7 20; 7 32; 7 36; 7 46; 8 4; 8 10; 8 15; 8 19; 8 26; 8 35; 8 47; 9 6; 9 11; 9 16; 9 39; 9 43; 9 45; 9 50; 10 9; 10 21; 10 32; 10 39; 11 4; 11 14; 11 15; 11 23; 11 32; 11 50; 12 21; 12 23; 13 8; 13 35; 13 38; 13 40; 13 49; 14 11; 14 16; 14 20; 14 31; 14 50; 15 13; 15 25; 15 37; 15 40; 16 2; 16 8; 16 21; 16 26; 16 32; 16 38; 16 41; 16 48; 17 6; 17 20; 17 22; 17 50; 18 5; 18 11; 18 21; 18 22; 18 26; 18 31; 18 49; 19 8; 19 21; 20 5; 20 14; 20 27; 20 39; 20 41; 20 44; 21 2; 21 4; 21 24; 21 37; 21 43; 22 9; 23 4; 23 8; 23 10; 23 11; 23 32; 23 40; 24 13; 24 14; 24 16; 24 27; 24 28; 24 34; 24 36; 25 4; 25 5; 25 8; 25 21; 25 45; 26 20; 26 23; 27 4; 27 48; 28 24; 28 39; 29 9; 29 18; 29 36; 30 15; 30 18; 30 25; 30 38; 30 46; 31 5; 31 6; 31 7; 31 28; 31 48; 32 6; 32 10; 32 27; 32 40; 32 49; 33 5; 33 13; 33 18; 33 22; 33 23; 33 26; 33 43; 33 46; 34 5; 34 9; 34 17; 34 27; 35 23; 35 24; 35 39; 35 41; 36 13; 36 40; 36 45; 37 8; 37 19; 37 28; 37 48; 38 11; 38 16; 38 17; 38 20; 38 32; 39 7; 39 8; 39 12; 39 28; 39 38; 39 44; 39 46; 39 50; 40 6; 40 14; 40 30; 40 45; 41 4; 41 15; 41 25; 41 27; 41 34; 42 10; 42 17; 42 24; 43 6; 43 44; 44 8; 44 21; 44 33; 44 39; 45 5; 45 10; 45 50; 46 5; 46 6; 46 7; 46 22; 46 27; 46 50; 47 22; 47 23; 47 24; 47 29; 47 39; 48 9; 48 14; 48 18; 48 36; 49 19; 49 24; 49 26; 49 27]
cL_orig = [21.0, 54.0, 71.0, 152.0, 246.0, 292.0, 316.0, 41.0, 71.0, 89.0, 155.0, 335.0, 437.0, 154.0, 185.0, 240.0, 255.0, 363.0, 413.0, 230.0, 276.0, 346.0, 363.0, 28.0, 149.0, 358.0, 366.0, 436.0, 452.0, 121.0, 394.0, 395.0, 104.0, 133.0, 245.0, 287.0, 395.0, 45.0, 7.0, 71.0, 111.0, 176.0, 266.0, 387.0, 27.0, 23.0, 73.0, 293.0, 338.0, 356.0, 415.0, 15.0, 115.0, 224.0, 286.0, 68.0, 28.0, 42.0, 125.0, 211.0, 391.0, 89.0, 114.0, 52.0, 211.0, 247.0, 270.0, 355.0, 29.0, 21.0, 64.0, 172.0, 363.0, 11.0, 96.0, 223.0, 247.0, 138.0, 79.0, 36.0, 102.0, 165.0, 216.0, 237.0, 317.0, 111.0, 29.0, 48.0, 330.0, 129.0, 68.0, 35.0, 43.0, 76.0, 130.0, 306.0, 111.0, 19.0, 151.0, 65.0, 70.0, 190.0, 215.0, 239.0, 189.0, 174.0, 28.0, 156.0, 225.0, 119.0, 185.0, 153.0, 125.0, 120.0, 88.0, 165.0, 112.0, 97.0, 70.0, 18.0, 38.0, 97.0, 112.0, 205.0, 190.0, 168.0, 30.0, 198.0, 61.0, 27.0, 233.0, 208.0, 42.0, 108.0, 198.0, 110.0, 59.0, 150.0, 118.0, 55.0, 78.0, 155.0, 253.0, 246.0, 236.0, 32.0, 166.0, 264.0, 219.0, 53.0, 78.0, 171.0, 280.0, 195.0, 152.0, 114.0, 104.0, 58.0, 96.0, 134.0, 292.0, 254.0, 168.0, 67.0, 124.0, 111.0, 43.0, 54.0, 234.0, 42.0, 81.0, 292.0, 166.0, 93.0, 104.0, 275.0, 218.0, 215.0, 176.0, 57.0, 315.0, 309.0, 270.0, 111.0, 7.0, 43.0, 70.0, 110.0, 330.0, 262.0, 99.0, 41.0, 359.0, 257.0, 163.0, 138.0, 71.0, 316.0, 249.0, 172.0, 365.0, 14.0, 361.0, 229.0, 98.0, 50.0, 400.0, 348.0, 52.0, 405.0, 399.0, 385.0, 239.0, 187.0, 30.0, 244.0, 239.0, 235.0, 184.0, 76.0, 390.0, 339.0, 294.0, 125.0, 304.0, 250.0, 227.0, 216.0]
cU_orig = [21.0, 64.0, 71.0, 152.0, 246.0, 310.0, 316.0, 41.0, 71.0, 89.0, 155.0, 335.0, 437.0, 168.0, 185.0, 240.0, 255.0, 363.0, 413.0, 230.0, 276.0, 356.0, 363.0, 28.0, 149.0, 358.0, 366.0, 436.0, 452.0, 121.0, 394.0, 395.0, 116.0, 133.0, 245.0, 303.0, 395.0, 45.0, 25.0, 71.0, 111.0, 176.0, 266.0, 387.0, 27.0, 23.0, 73.0, 309.0, 338.0, 356.0, 415.0, 15.0, 115.0, 224.0, 286.0, 82.0, 38.0, 42.0, 125.0, 211.0, 391.0, 89.0, 114.0, 52.0, 221.0, 247.0, 270.0, 355.0, 29.0, 21.0, 64.0, 172.0, 363.0, 29.0, 96.0, 223.0, 247.0, 138.0, 79.0, 56.0, 102.0, 165.0, 230.0, 255.0, 333.0, 111.0, 29.0, 48.0, 330.0, 129.0, 68.0, 35.0, 43.0, 76.0, 130.0, 306.0, 111.0, 19.0, 151.0, 65.0, 70.0, 190.0, 215.0, 239.0, 189.0, 174.0, 28.0, 156.0, 225.0, 133.0, 185.0, 153.0, 125.0, 120.0, 88.0, 177.0, 112.0, 97.0, 88.0, 36.0, 38.0, 97.0, 124.0, 205.0, 202.0, 168.0, 42.0, 198.0, 61.0, 27.0, 233.0, 208.0, 42.0, 108.0, 198.0, 110.0, 79.0, 150.0, 118.0, 55.0, 78.0, 155.0, 263.0, 264.0, 236.0, 32.0, 166.0, 264.0, 219.0, 53.0, 90.0, 171.0, 280.0, 195.0, 152.0, 114.0, 104.0, 74.0, 96.0, 134.0, 292.0, 254.0, 168.0, 67.0, 124.0, 111.0, 43.0, 74.0, 234.0, 42.0, 99.0, 292.0, 186.0, 93.0, 120.0, 275.0, 218.0, 215.0, 176.0, 73.0, 315.0, 309.0, 270.0, 111.0, 21.0, 59.0, 70.0, 110.0, 342.0, 262.0, 99.0, 53.0, 371.0, 257.0, 163.0, 138.0, 71.0, 316.0, 249.0, 188.0, 365.0, 14.0, 361.0, 229.0, 112.0, 50.0, 400.0, 348.0, 52.0, 405.0, 411.0, 385.0, 239.0, 187.0, 42.0, 256.0, 251.0, 235.0, 184.0, 76.0, 390.0, 339.0, 310.0, 125.0, 304.0, 250.0, 227.0, 232.0]
d = [14.0, 6.0, 16.0, 16.0, 15.0, 20.0, 16.0, 10.0, 7.0, 11.0, 2.0, 12.0, 7.0, 2.0, 5.0, 15.0, 14.0, 4.0, 12.0, 16.0, 4.0, 16.0, 5.0, 20.0, 1.0, 2.0, 4.0, 5.0, 17.0, 7.0, 15.0, 5.0, 6.0, 19.0, 6.0, 18.0, 19.0, 8.0, 9.0, 1.0, 5.0, 14.0, 10.0, 20.0, 5.0, 18.0, 14.0, 19.0, 17.0, 6.0, 16.0, 15.0, 9.0, 8.0, 16.0, 6.0, 14.0, 1.0, 2.0, 12.0, 9.0, 5.0, 3.0, 12.0, 11.0, 6.0, 19.0, 1.0, 13.0, 10.0, 7.0, 8.0, 13.0, 12.0, 8.0, 6.0, 15.0, 7.0, 2.0, 17.0, 6.0, 8.0, 4.0, 8.0, 13.0, 19.0, 12.0, 9.0, 5.0, 16.0, 7.0, 6.0, 12.0, 15.0, 16.0, 6.0, 4.0, 17.0, 19.0, 17.0, 2.0, 4.0, 1.0, 15.0, 4.0, 19.0, 19.0, 19.0, 20.0, 20.0, 19.0, 20.0, 11.0, 4.0, 15.0, 4.0, 12.0, 16.0, 15.0, 15.0, 19.0, 14.0, 19.0, 8.0, 6.0, 13.0, 15.0, 6.0, 12.0, 9.0, 7.0, 8.0, 9.0, 19.0, 7.0, 16.0, 19.0, 19.0, 14.0, 16.0, 11.0, 18.0, 17.0, 5.0, 19.0, 12.0, 6.0, 6.0, 4.0, 19.0, 8.0, 12.0, 18.0, 18.0, 4.0, 12.0, 2.0, 1.0, 17.0, 9.0, 6.0, 3.0, 19.0, 3.0, 15.0, 10.0, 17.0, 15.0, 12.0, 20.0, 17.0, 17.0, 10.0, 20.0, 4.0, 11.0, 16.0, 6.0, 19.0, 7.0, 16.0, 9.0, 18.0, 14.0, 3.0, 13.0, 9.0, 20.0, 10.0, 11.0, 14.0, 20.0, 20.0, 9.0, 9.0, 14.0, 7.0, 14.0, 14.0, 11.0, 15.0, 2.0, 20.0, 19.0, 6.0, 5.0, 1.0, 20.0, 19.0, 8.0, 15.0, 8.0, 4.0, 13.0, 20.0, 8.0, 13.0, 6.0, 5.0, 8.0, 17.0, 1.0, 8.0, 14.0, 11.0, 8.0, 4.0, 19.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.023
