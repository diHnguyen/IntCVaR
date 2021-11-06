edge = [1 13; 1 43; 1 46; 2 4; 2 7; 2 12; 2 33; 2 44; 3 10; 3 12; 3 29; 3 41; 4 19; 4 43; 4 45; 5 4; 5 8; 5 18; 5 37; 6 7; 6 8; 6 23; 6 26; 6 35; 6 44; 7 4; 7 12; 7 25; 7 31; 7 44; 8 2; 8 11; 8 27; 9 4; 9 5; 9 6; 9 13; 9 32; 9 45; 10 6; 10 14; 10 44; 11 23; 11 35; 11 37; 12 15; 12 16; 12 20; 12 38; 13 7; 13 8; 13 18; 13 21; 13 30; 14 6; 14 21; 14 23; 14 29; 14 36; 14 47; 15 4; 15 20; 15 22; 15 24; 15 26; 15 30; 15 31; 15 32; 15 41; 16 3; 16 15; 16 47; 17 8; 17 36; 17 41; 17 42; 17 47; 18 8; 18 15; 18 21; 18 26; 18 28; 18 31; 18 47; 19 9; 19 25; 19 29; 19 33; 19 36; 20 18; 20 19; 20 35; 20 39; 21 9; 21 12; 21 15; 21 16; 21 27; 21 38; 22 8; 22 15; 22 17; 22 32; 22 34; 22 35; 23 16; 23 25; 23 27; 23 33; 23 41; 24 8; 24 23; 24 36; 25 17; 25 26; 26 10; 26 11; 26 33; 26 49; 27 46; 27 50; 28 37; 28 42; 28 47; 29 4; 29 6; 29 13; 29 18; 29 38; 29 45; 30 7; 30 21; 30 25; 30 41; 30 42; 30 43; 30 45; 30 49; 31 5; 31 8; 31 35; 31 39; 31 50; 32 12; 32 38; 32 44; 32 47; 33 14; 33 30; 33 36; 33 50; 34 14; 34 19; 34 28; 34 33; 34 48; 34 49; 35 6; 35 11; 35 25; 35 27; 35 28; 35 36; 35 42; 36 33; 36 40; 36 47; 37 2; 37 4; 37 14; 37 22; 37 29; 37 31; 37 38; 37 39; 38 17; 38 32; 38 35; 38 37; 39 4; 39 12; 39 33; 39 50; 40 15; 40 25; 40 26; 40 45; 41 8; 41 13; 41 32; 41 36; 41 44; 42 13; 42 15; 42 39; 42 43; 43 6; 43 42; 44 10; 44 12; 44 15; 44 32; 44 41; 44 48; 45 3; 45 24; 45 47; 46 3; 46 17; 46 28; 46 45; 47 20; 48 8; 48 11; 48 22; 48 33; 48 39; 48 46; 49 11; 49 18; 49 30; 49 31]
cL_orig = [114.0, 421.0, 447.0, 17.0, 55.0, 97.0, 311.0, 413.0, 67.0, 89.0, 265.0, 382.0, 146.0, 389.0, 405.0, 13.0, 29.0, 125.0, 324.0, 6.0, 19.0, 175.0, 187.0, 276.0, 377.0, 34.0, 55.0, 181.0, 242.0, 371.0, 55.0, 28.0, 179.0, 52.0, 36.0, 25.0, 35.0, 224.0, 364.0, 37.0, 36.0, 335.0, 121.0, 243.0, 255.0, 25.0, 39.0, 79.0, 257.0, 63.0, 48.0, 46.0, 79.0, 173.0, 76.0, 61.0, 85.0, 147.0, 222.0, 327.0, 100.0, 37.0, 69.0, 89.0, 107.0, 151.0, 155.0, 168.0, 257.0, 127.0, 6.0, 315.0, 76.0, 191.0, 230.0, 240.0, 301.0, 95.0, 30.0, 32.0, 85.0, 96.0, 117.0, 288.0, 104.0, 62.0, 99.0, 144.0, 173.0, 15.0, 11.0, 139.0, 189.0, 109.0, 82.0, 59.0, 46.0, 65.0, 160.0, 144.0, 65.0, 49.0, 91.0, 119.0, 128.0, 63.0, 19.0, 36.0, 100.0, 178.0, 164.0, 5.0, 121.0, 65.0, 0.0, 155.0, 149.0, 67.0, 216.0, 185.0, 228.0, 86.0, 143.0, 185.0, 249.0, 228.0, 154.0, 114.0, 88.0, 162.0, 226.0, 90.0, 52.0, 106.0, 118.0, 130.0, 145.0, 180.0, 256.0, 230.0, 38.0, 79.0, 178.0, 204.0, 59.0, 116.0, 136.0, 194.0, 30.0, 18.0, 169.0, 195.0, 146.0, 60.0, 0.0, 141.0, 154.0, 287.0, 236.0, 98.0, 78.0, 67.0, 0.0, 65.0, 35.0, 42.0, 111.0, 347.0, 332.0, 232.0, 148.0, 80.0, 57.0, 5.0, 19.0, 206.0, 58.0, 23.0, 5.0, 352.0, 266.0, 62.0, 108.0, 247.0, 142.0, 140.0, 55.0, 329.0, 270.0, 95.0, 47.0, 31.0, 295.0, 274.0, 25.0, 14.0, 375.0, 14.0, 338.0, 321.0, 292.0, 113.0, 31.0, 36.0, 420.0, 215.0, 15.0, 428.0, 294.0, 182.0, 8.0, 274.0, 403.0, 366.0, 262.0, 141.0, 88.0, 20.0, 381.0, 309.0, 189.0, 179.0]
cU_orig = [134.0, 421.0, 447.0, 17.0, 55.0, 97.0, 311.0, 425.0, 67.0, 89.0, 265.0, 382.0, 146.0, 389.0, 405.0, 13.0, 29.0, 125.0, 324.0, 6.0, 19.0, 175.0, 207.0, 294.0, 389.0, 34.0, 55.0, 181.0, 242.0, 371.0, 55.0, 28.0, 193.0, 52.0, 36.0, 25.0, 35.0, 244.0, 364.0, 37.0, 36.0, 335.0, 121.0, 243.0, 255.0, 43.0, 39.0, 79.0, 257.0, 63.0, 48.0, 56.0, 79.0, 173.0, 76.0, 77.0, 85.0, 147.0, 222.0, 327.0, 114.0, 57.0, 69.0, 89.0, 121.0, 151.0, 167.0, 168.0, 257.0, 127.0, 6.0, 315.0, 96.0, 191.0, 248.0, 260.0, 301.0, 95.0, 40.0, 32.0, 85.0, 96.0, 137.0, 288.0, 104.0, 62.0, 99.0, 144.0, 173.0, 15.0, 11.0, 151.0, 189.0, 127.0, 98.0, 59.0, 56.0, 65.0, 180.0, 144.0, 65.0, 61.0, 105.0, 119.0, 128.0, 75.0, 19.0, 36.0, 100.0, 178.0, 164.0, 5.0, 121.0, 85.0, 18.0, 155.0, 149.0, 67.0, 234.0, 185.0, 228.0, 86.0, 143.0, 185.0, 249.0, 228.0, 174.0, 114.0, 88.0, 162.0, 226.0, 90.0, 52.0, 106.0, 118.0, 130.0, 145.0, 190.0, 256.0, 230.0, 38.0, 91.0, 194.0, 204.0, 59.0, 116.0, 156.0, 194.0, 30.0, 34.0, 169.0, 195.0, 146.0, 60.0, 20.0, 141.0, 154.0, 287.0, 236.0, 98.0, 78.0, 67.0, 20.0, 65.0, 35.0, 42.0, 111.0, 363.0, 332.0, 232.0, 148.0, 80.0, 57.0, 5.0, 19.0, 206.0, 58.0, 37.0, 5.0, 352.0, 266.0, 62.0, 108.0, 247.0, 162.0, 140.0, 55.0, 329.0, 280.0, 95.0, 47.0, 31.0, 295.0, 274.0, 25.0, 14.0, 375.0, 14.0, 338.0, 321.0, 292.0, 123.0, 31.0, 36.0, 420.0, 215.0, 15.0, 428.0, 294.0, 182.0, 8.0, 274.0, 403.0, 366.0, 262.0, 153.0, 88.0, 20.0, 381.0, 309.0, 189.0, 179.0]
d = [7.0, 2.0, 11.0, 12.0, 19.0, 18.0, 17.0, 9.0, 8.0, 19.0, 11.0, 13.0, 3.0, 11.0, 3.0, 13.0, 9.0, 6.0, 1.0, 5.0, 2.0, 14.0, 15.0, 12.0, 18.0, 2.0, 13.0, 14.0, 14.0, 5.0, 15.0, 8.0, 10.0, 3.0, 19.0, 12.0, 1.0, 15.0, 11.0, 11.0, 17.0, 8.0, 13.0, 20.0, 20.0, 5.0, 12.0, 10.0, 3.0, 7.0, 19.0, 16.0, 7.0, 2.0, 17.0, 11.0, 11.0, 3.0, 7.0, 5.0, 14.0, 17.0, 6.0, 13.0, 5.0, 9.0, 13.0, 17.0, 18.0, 12.0, 9.0, 5.0, 1.0, 9.0, 9.0, 9.0, 7.0, 17.0, 16.0, 4.0, 7.0, 20.0, 6.0, 15.0, 20.0, 19.0, 7.0, 1.0, 14.0, 7.0, 19.0, 19.0, 17.0, 1.0, 14.0, 14.0, 17.0, 20.0, 17.0, 13.0, 15.0, 18.0, 17.0, 5.0, 10.0, 4.0, 12.0, 15.0, 13.0, 13.0, 17.0, 11.0, 13.0, 11.0, 19.0, 2.0, 7.0, 17.0, 7.0, 16.0, 20.0, 5.0, 5.0, 10.0, 7.0, 5.0, 14.0, 20.0, 18.0, 11.0, 3.0, 16.0, 13.0, 11.0, 7.0, 17.0, 11.0, 13.0, 6.0, 8.0, 1.0, 13.0, 9.0, 15.0, 18.0, 10.0, 18.0, 7.0, 6.0, 15.0, 2.0, 11.0, 15.0, 7.0, 16.0, 17.0, 9.0, 20.0, 6.0, 6.0, 5.0, 3.0, 1.0, 18.0, 3.0, 14.0, 19.0, 2.0, 3.0, 14.0, 12.0, 18.0, 20.0, 12.0, 5.0, 1.0, 14.0, 14.0, 13.0, 20.0, 6.0, 10.0, 20.0, 6.0, 4.0, 20.0, 10.0, 19.0, 3.0, 8.0, 3.0, 14.0, 9.0, 9.0, 9.0, 17.0, 6.0, 20.0, 16.0, 4.0, 20.0, 10.0, 15.0, 12.0, 5.0, 6.0, 4.0, 18.0, 4.0, 3.0, 14.0, 1.0, 10.0, 14.0, 8.0, 9.0, 8.0, 16.0, 11.0, 8.0, 15.0, 3.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
b = 4
last_node = maximum(edge)
β = 0.821
