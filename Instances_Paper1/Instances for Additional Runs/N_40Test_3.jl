edge = [1 3; 1 6; 1 8; 1 12; 1 15; 1 30; 1 33; 2 5; 2 9; 2 10; 2 25; 2 27; 2 29; 2 35; 2 40; 3 5; 3 9; 3 28; 3 31; 3 35; 3 37; 3 40; 4 2; 4 5; 4 9; 4 11; 4 12; 4 15; 4 16; 4 19; 4 22; 4 23; 4 25; 4 26; 4 31; 4 32; 4 35; 4 36; 4 37; 5 15; 5 29; 5 34; 5 40; 6 2; 6 3; 6 9; 6 18; 6 21; 6 24; 6 36; 6 37; 7 6; 7 13; 7 22; 7 26; 7 28; 7 30; 7 34; 7 40; 8 16; 8 18; 8 20; 8 24; 8 38; 9 2; 9 10; 9 11; 9 13; 9 24; 10 3; 10 4; 10 6; 10 9; 10 13; 10 18; 10 22; 10 27; 10 28; 10 38; 10 39; 10 40; 11 5; 11 6; 11 21; 11 33; 11 35; 11 39; 11 40; 12 9; 12 19; 12 24; 12 26; 12 30; 12 34; 12 37; 13 16; 13 24; 13 27; 14 3; 14 17; 14 28; 14 30; 14 37; 14 40; 15 5; 15 9; 15 12; 15 14; 15 17; 15 21; 15 36; 15 37; 15 38; 15 39; 16 11; 16 14; 16 15; 16 19; 16 20; 16 34; 16 37; 17 4; 17 11; 17 15; 17 21; 17 23; 17 24; 17 25; 17 27; 17 39; 17 40; 18 4; 18 6; 18 13; 18 16; 18 27; 18 35; 18 38; 19 2; 19 3; 19 5; 19 6; 19 9; 19 11; 19 25; 19 26; 19 37; 19 38; 20 5; 20 12; 20 15; 20 22; 20 26; 20 37; 21 5; 21 17; 21 25; 21 28; 21 30; 21 31; 22 4; 22 10; 22 23; 22 29; 22 30; 22 32; 22 38; 22 40; 23 4; 23 9; 23 24; 23 28; 23 29; 23 36; 23 39; 24 3; 24 5; 24 21; 24 22; 24 29; 24 31; 24 40; 25 2; 25 7; 25 10; 25 14; 25 19; 25 23; 25 24; 25 27; 25 31; 25 32; 26 14; 26 17; 26 24; 26 27; 26 33; 26 36; 26 38; 26 39; 27 6; 27 9; 27 14; 27 20; 27 24; 27 29; 27 35; 27 36; 27 37; 28 6; 28 8; 28 9; 28 16; 28 38; 29 2; 29 4; 29 10; 29 15; 29 18; 29 21; 29 23; 29 24; 29 37; 29 40; 30 6; 30 9; 30 10; 30 14; 30 24; 30 26; 30 38; 30 39; 31 3; 31 5; 31 6; 31 16; 31 18; 31 19; 31 20; 31 29; 31 32; 31 36; 32 2; 32 4; 32 5; 32 13; 32 14; 32 16; 32 19; 32 37; 32 38; 33 9; 33 10; 33 14; 33 15; 33 32; 33 38; 33 39; 34 4; 34 22; 34 26; 34 29; 34 36; 34 37; 34 38; 35 3; 35 4; 35 6; 35 12; 35 22; 35 25; 35 31; 35 39; 36 3; 36 8; 36 10; 36 12; 36 13; 36 16; 36 32; 37 2; 37 4; 37 15; 37 29; 37 38; 38 2; 38 5; 38 6; 38 10; 38 17; 38 34; 38 36; 38 40; 39 2; 39 3; 39 9; 39 11; 39 18; 39 19; 39 22; 39 27; 39 30; 39 33; 39 36]
cL_orig = [19.0, 47.0, 65.0, 114.0, 139.0, 293.0, 319.0, 34.0, 59.0, 72.0, 228.0, 247.0, 272.0, 333.0, 382.0, 19.0, 58.0, 250.0, 281.0, 319.0, 342.0, 372.0, 25.0, 11.0, 50.0, 74.0, 85.0, 115.0, 119.0, 145.0, 178.0, 195.0, 213.0, 217.0, 270.0, 278.0, 309.0, 324.0, 332.0, 95.0, 241.0, 292.0, 349.0, 42.0, 32.0, 27.0, 120.0, 138.0, 177.0, 297.0, 310.0, 9.0, 61.0, 149.0, 187.0, 208.0, 233.0, 274.0, 330.0, 78.0, 96.0, 120.0, 159.0, 305.0, 75.0, 10.0, 21.0, 41.0, 149.0, 72.0, 58.0, 40.0, 12.0, 32.0, 84.0, 121.0, 167.0, 182.0, 279.0, 294.0, 299.0, 58.0, 48.0, 98.0, 219.0, 243.0, 283.0, 285.0, 28.0, 67.0, 120.0, 145.0, 179.0, 215.0, 252.0, 31.0, 106.0, 144.0, 108.0, 25.0, 143.0, 156.0, 233.0, 257.0, 82.0, 64.0, 32.0, 5.0, 19.0, 61.0, 213.0, 225.0, 232.0, 244.0, 45.0, 25.0, 11.0, 27.0, 39.0, 183.0, 208.0, 129.0, 55.0, 15.0, 37.0, 58.0, 71.0, 81.0, 100.0, 223.0, 234.0, 144.0, 125.0, 47.0, 21.0, 90.0, 168.0, 198.0, 166.0, 163.0, 140.0, 134.0, 95.0, 80.0, 63.0, 71.0, 181.0, 185.0, 148.0, 82.0, 50.0, 20.0, 56.0, 173.0, 157.0, 42.0, 43.0, 65.0, 91.0, 99.0, 180.0, 118.0, 6.0, 66.0, 78.0, 104.0, 159.0, 181.0, 190.0, 144.0, 12.0, 54.0, 63.0, 126.0, 164.0, 212.0, 190.0, 29.0, 18.0, 47.0, 73.0, 163.0, 225.0, 185.0, 145.0, 113.0, 63.0, 12.0, 15.0, 10.0, 56.0, 67.0, 122.0, 91.0, 16.0, 6.0, 72.0, 96.0, 117.0, 133.0, 210.0, 182.0, 128.0, 75.0, 31.0, 21.0, 79.0, 91.0, 102.0, 217.0, 201.0, 191.0, 119.0, 97.0, 271.0, 247.0, 191.0, 139.0, 110.0, 84.0, 62.0, 47.0, 85.0, 112.0, 235.0, 205.0, 195.0, 162.0, 56.0, 31.0, 78.0, 87.0, 279.0, 265.0, 253.0, 145.0, 131.0, 117.0, 112.0, 22.0, 7.0, 45.0, 304.0, 283.0, 268.0, 188.0, 177.0, 164.0, 128.0, 51.0, 56.0, 245.0, 230.0, 193.0, 181.0, 10.0, 52.0, 58.0, 301.0, 123.0, 85.0, 51.0, 18.0, 30.0, 42.0, 315.0, 313.0, 291.0, 226.0, 126.0, 81.0, 41.0, 27.0, 333.0, 275.0, 259.0, 244.0, 228.0, 196.0, 40.0, 345.0, 328.0, 220.0, 77.0, 6.0, 356.0, 332.0, 323.0, 276.0, 210.0, 42.0, 20.0, 21.0, 371.0, 364.0, 295.0, 285.0, 210.0, 196.0, 172.0, 124.0, 92.0, 65.0, 27.0]
cU_orig = [19.0, 47.0, 65.0, 114.0, 139.0, 293.0, 319.0, 34.0, 79.0, 88.0, 228.0, 247.0, 272.0, 333.0, 382.0, 19.0, 58.0, 250.0, 281.0, 319.0, 342.0, 372.0, 25.0, 11.0, 50.0, 74.0, 85.0, 115.0, 119.0, 145.0, 178.0, 195.0, 213.0, 217.0, 270.0, 278.0, 309.0, 324.0, 332.0, 95.0, 241.0, 292.0, 349.0, 42.0, 32.0, 27.0, 120.0, 164.0, 177.0, 297.0, 310.0, 9.0, 61.0, 149.0, 187.0, 208.0, 233.0, 274.0, 330.0, 78.0, 96.0, 120.0, 159.0, 305.0, 75.0, 10.0, 21.0, 41.0, 149.0, 72.0, 58.0, 40.0, 12.0, 32.0, 84.0, 121.0, 167.0, 182.0, 279.0, 294.0, 299.0, 58.0, 48.0, 98.0, 219.0, 243.0, 283.0, 285.0, 28.0, 67.0, 120.0, 145.0, 179.0, 215.0, 252.0, 31.0, 106.0, 144.0, 108.0, 25.0, 143.0, 156.0, 233.0, 257.0, 112.0, 64.0, 32.0, 5.0, 19.0, 61.0, 213.0, 225.0, 232.0, 244.0, 45.0, 25.0, 11.0, 27.0, 39.0, 183.0, 208.0, 129.0, 55.0, 15.0, 37.0, 58.0, 71.0, 81.0, 100.0, 223.0, 234.0, 144.0, 125.0, 47.0, 21.0, 90.0, 178.0, 198.0, 166.0, 163.0, 140.0, 134.0, 95.0, 80.0, 63.0, 71.0, 181.0, 185.0, 148.0, 82.0, 50.0, 20.0, 56.0, 173.0, 157.0, 42.0, 43.0, 65.0, 91.0, 99.0, 180.0, 118.0, 6.0, 66.0, 78.0, 104.0, 159.0, 181.0, 190.0, 144.0, 12.0, 54.0, 63.0, 126.0, 164.0, 212.0, 190.0, 29.0, 18.0, 47.0, 73.0, 163.0, 225.0, 185.0, 145.0, 113.0, 63.0, 38.0, 15.0, 28.0, 56.0, 67.0, 122.0, 91.0, 16.0, 6.0, 72.0, 96.0, 117.0, 133.0, 210.0, 182.0, 128.0, 75.0, 31.0, 21.0, 79.0, 91.0, 102.0, 217.0, 201.0, 191.0, 119.0, 97.0, 271.0, 247.0, 191.0, 139.0, 110.0, 84.0, 62.0, 47.0, 85.0, 112.0, 235.0, 205.0, 195.0, 162.0, 56.0, 45.0, 78.0, 87.0, 279.0, 265.0, 253.0, 145.0, 131.0, 117.0, 112.0, 22.0, 7.0, 45.0, 304.0, 283.0, 268.0, 188.0, 177.0, 164.0, 128.0, 51.0, 56.0, 245.0, 230.0, 193.0, 181.0, 10.0, 52.0, 58.0, 301.0, 123.0, 85.0, 51.0, 18.0, 30.0, 42.0, 315.0, 313.0, 291.0, 226.0, 126.0, 109.0, 41.0, 47.0, 333.0, 275.0, 259.0, 244.0, 228.0, 196.0, 40.0, 345.0, 328.0, 220.0, 77.0, 6.0, 356.0, 332.0, 323.0, 276.0, 210.0, 42.0, 20.0, 21.0, 371.0, 364.0, 295.0, 285.0, 210.0, 196.0, 172.0, 124.0, 92.0, 65.0, 27.0]
d = [2.0, 5.0, 7.0, 12.0, 17.0, 19.0, 3.0, 18.0, 4.0, 11.0, 10.0, 8.0, 9.0, 2.0, 13.0, 12.0, 17.0, 3.0, 5.0, 6.0, 14.0, 13.0, 20.0, 1.0, 20.0, 20.0, 16.0, 14.0, 15.0, 1.0, 6.0, 9.0, 16.0, 17.0, 14.0, 4.0, 6.0, 20.0, 1.0, 16.0, 15.0, 10.0, 5.0, 7.0, 14.0, 1.0, 18.0, 1.0, 7.0, 11.0, 7.0, 17.0, 12.0, 6.0, 19.0, 5.0, 10.0, 13.0, 13.0, 19.0, 18.0, 5.0, 8.0, 2.0, 14.0, 17.0, 7.0, 17.0, 14.0, 14.0, 11.0, 6.0, 14.0, 6.0, 19.0, 13.0, 3.0, 5.0, 15.0, 2.0, 15.0, 4.0, 6.0, 20.0, 13.0, 8.0, 2.0, 4.0, 14.0, 9.0, 6.0, 12.0, 2.0, 6.0, 4.0, 18.0, 18.0, 6.0, 1.0, 1.0, 1.0, 10.0, 11.0, 18.0, 16.0, 9.0, 7.0, 20.0, 1.0, 9.0, 8.0, 18.0, 8.0, 8.0, 1.0, 14.0, 16.0, 20.0, 16.0, 5.0, 18.0, 11.0, 1.0, 14.0, 19.0, 20.0, 10.0, 2.0, 2.0, 20.0, 3.0, 5.0, 17.0, 9.0, 16.0, 5.0, 19.0, 18.0, 9.0, 2.0, 17.0, 16.0, 7.0, 17.0, 20.0, 15.0, 12.0, 2.0, 10.0, 13.0, 9.0, 14.0, 8.0, 15.0, 17.0, 6.0, 20.0, 7.0, 13.0, 15.0, 3.0, 16.0, 12.0, 11.0, 1.0, 17.0, 12.0, 14.0, 2.0, 4.0, 17.0, 3.0, 18.0, 9.0, 19.0, 11.0, 16.0, 6.0, 6.0, 14.0, 1.0, 20.0, 12.0, 1.0, 11.0, 3.0, 4.0, 7.0, 16.0, 13.0, 16.0, 2.0, 18.0, 15.0, 11.0, 20.0, 11.0, 5.0, 15.0, 17.0, 15.0, 17.0, 20.0, 9.0, 8.0, 14.0, 9.0, 5.0, 20.0, 7.0, 17.0, 14.0, 6.0, 3.0, 11.0, 2.0, 18.0, 12.0, 7.0, 1.0, 15.0, 10.0, 8.0, 12.0, 10.0, 12.0, 14.0, 13.0, 10.0, 11.0, 8.0, 14.0, 17.0, 13.0, 1.0, 13.0, 14.0, 17.0, 3.0, 18.0, 6.0, 9.0, 8.0, 4.0, 1.0, 9.0, 4.0, 17.0, 5.0, 8.0, 20.0, 16.0, 17.0, 2.0, 13.0, 9.0, 7.0, 10.0, 13.0, 2.0, 12.0, 5.0, 2.0, 1.0, 14.0, 4.0, 10.0, 2.0, 11.0, 6.0, 10.0, 12.0, 15.0, 2.0, 20.0, 4.0, 13.0, 2.0, 15.0, 8.0, 4.0, 13.0, 4.0, 1.0, 16.0, 19.0, 10.0, 18.0, 16.0, 17.0, 1.0, 4.0, 20.0, 20.0, 7.0, 16.0, 18.0, 11.0, 13.0, 7.0, 19.0, 14.0, 18.0, 20.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 40
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 5.0
b = 2
Grp = center
last_node = maximum(edge)