edge = [1 4; 1 12; 1 23; 1 36; 1 46; 2 3; 2 6; 2 16; 2 19; 2 21; 2 36; 3 33; 3 46; 4 5; 4 8; 4 9; 4 16; 4 19; 4 22; 4 27; 4 28; 4 30; 4 31; 5 2; 5 21; 5 22; 5 26; 5 27; 6 13; 6 18; 6 26; 7 15; 7 27; 7 28; 7 44; 8 6; 8 14; 8 23; 8 28; 8 31; 8 40; 8 46; 8 49; 9 26; 9 35; 9 44; 9 49; 9 50; 10 7; 10 20; 10 23; 10 37; 11 16; 11 17; 11 29; 11 37; 12 23; 12 33; 12 40; 13 17; 13 32; 13 41; 14 6; 14 7; 14 23; 14 35; 14 38; 14 43; 15 11; 15 39; 16 2; 16 24; 16 36; 16 38; 16 42; 17 31; 17 39; 18 13; 18 16; 18 20; 18 25; 18 27; 18 34; 18 41; 18 49; 19 8; 19 9; 19 34; 19 40; 19 47; 20 7; 20 23; 20 39; 20 42; 21 5; 21 20; 21 26; 21 35; 21 38; 21 43; 21 47; 22 7; 22 21; 22 32; 23 15; 23 41; 23 45; 23 50; 24 10; 24 14; 24 48; 24 50; 25 3; 25 5; 25 6; 25 9; 25 10; 25 11; 25 16; 26 15; 26 16; 26 25; 26 39; 26 47; 27 21; 27 39; 27 40; 27 46; 28 2; 28 9; 28 11; 28 12; 28 14; 28 17; 28 25; 28 26; 28 32; 29 18; 29 47; 29 49; 29 50; 30 36; 31 3; 31 13; 31 15; 31 23; 31 25; 31 27; 31 35; 31 37; 31 42; 31 48; 32 12; 32 27; 32 28; 32 38; 32 45; 32 48; 33 23; 33 29; 33 39; 34 10; 34 16; 34 23; 34 31; 35 2; 35 18; 35 21; 35 28; 35 43; 36 42; 36 46; 36 50; 37 30; 37 48; 37 49; 38 2; 38 9; 38 28; 38 29; 38 30; 38 32; 38 34; 39 40; 39 42; 40 2; 40 4; 40 20; 40 28; 40 30; 40 33; 40 43; 41 12; 41 28; 41 29; 41 30; 41 34; 41 38; 42 16; 42 20; 42 24; 42 25; 42 28; 42 40; 42 43; 42 46; 43 6; 43 9; 43 29; 44 3; 44 12; 44 16; 44 27; 44 31; 44 48; 45 3; 45 12; 45 13; 45 32; 45 36; 46 4; 46 12; 46 22; 46 36; 46 49; 47 20; 47 33; 47 46; 48 14; 48 16; 48 20; 48 24; 48 25; 48 30; 48 31; 48 36; 48 39; 48 43; 48 44; 48 49; 49 5; 49 6; 49 8; 49 17; 49 19; 49 21; 49 26; 49 27; 49 37]
cL_orig = [33.0, 112.0, 213.0, 351.0, 447.0, 7.0, 41.0, 135.0, 172.0, 191.0, 341.0, 304.0, 431.0, 9.0, 42.0, 51.0, 125.0, 147.0, 179.0, 231.0, 238.0, 261.0, 265.0, 32.0, 161.0, 169.0, 207.0, 216.0, 75.0, 117.0, 205.0, 75.0, 201.0, 206.0, 371.0, 22.0, 65.0, 153.0, 195.0, 224.0, 325.0, 376.0, 407.0, 174.0, 264.0, 345.0, 402.0, 407.0, 26.0, 105.0, 132.0, 275.0, 45.0, 56.0, 184.0, 262.0, 107.0, 209.0, 281.0, 28.0, 185.0, 284.0, 81.0, 70.0, 95.0, 208.0, 238.0, 285.0, 44.0, 235.0, 139.0, 77.0, 204.0, 215.0, 263.0, 144.0, 213.0, 52.0, 16.0, 14.0, 71.0, 86.0, 156.0, 234.0, 314.0, 108.0, 95.0, 149.0, 208.0, 283.0, 133.0, 31.0, 192.0, 206.0, 164.0, 13.0, 51.0, 135.0, 166.0, 220.0, 261.0, 145.0, 15.0, 97.0, 77.0, 177.0, 217.0, 271.0, 144.0, 100.0, 242.0, 255.0, 224.0, 197.0, 194.0, 160.0, 151.0, 141.0, 90.0, 112.0, 102.0, 13.0, 125.0, 212.0, 60.0, 117.0, 131.0, 186.0, 258.0, 190.0, 165.0, 157.0, 134.0, 109.0, 30.0, 25.0, 37.0, 115.0, 183.0, 199.0, 210.0, 59.0, 281.0, 182.0, 156.0, 81.0, 50.0, 44.0, 43.0, 61.0, 106.0, 170.0, 198.0, 50.0, 36.0, 65.0, 133.0, 156.0, 95.0, 40.0, 64.0, 235.0, 176.0, 106.0, 26.0, 331.0, 175.0, 144.0, 75.0, 82.0, 61.0, 104.0, 134.0, 68.0, 108.0, 118.0, 360.0, 290.0, 95.0, 89.0, 85.0, 45.0, 43.0, 6.0, 30.0, 377.0, 359.0, 193.0, 116.0, 96.0, 70.0, 27.0, 291.0, 125.0, 118.0, 111.0, 65.0, 28.0, 259.0, 222.0, 182.0, 167.0, 142.0, 23.0, 15.0, 33.0, 373.0, 342.0, 141.0, 414.0, 317.0, 271.0, 167.0, 129.0, 43.0, 417.0, 332.0, 325.0, 132.0, 93.0, 416.0, 331.0, 241.0, 104.0, 33.0, 273.0, 144.0, 5.0, 338.0, 316.0, 276.0, 245.0, 220.0, 181.0, 166.0, 113.0, 93.0, 52.0, 35.0, 7.0, 435.0, 428.0, 406.0, 315.0, 296.0, 281.0, 228.0, 213.0, 122.0]
cU_orig = [33.0, 112.0, 225.0, 351.0, 447.0, 7.0, 41.0, 135.0, 172.0, 191.0, 341.0, 304.0, 431.0, 9.0, 42.0, 51.0, 125.0, 147.0, 179.0, 231.0, 238.0, 261.0, 265.0, 32.0, 161.0, 169.0, 207.0, 216.0, 75.0, 117.0, 205.0, 75.0, 201.0, 206.0, 371.0, 22.0, 65.0, 153.0, 195.0, 234.0, 325.0, 376.0, 407.0, 174.0, 264.0, 345.0, 402.0, 423.0, 40.0, 105.0, 132.0, 275.0, 55.0, 56.0, 184.0, 262.0, 107.0, 209.0, 281.0, 44.0, 185.0, 284.0, 81.0, 70.0, 95.0, 208.0, 250.0, 285.0, 44.0, 235.0, 139.0, 77.0, 204.0, 215.0, 263.0, 144.0, 233.0, 52.0, 16.0, 30.0, 71.0, 86.0, 156.0, 234.0, 314.0, 108.0, 95.0, 149.0, 208.0, 283.0, 133.0, 31.0, 192.0, 226.0, 164.0, 13.0, 51.0, 135.0, 166.0, 220.0, 261.0, 145.0, 15.0, 97.0, 77.0, 177.0, 217.0, 271.0, 144.0, 100.0, 242.0, 255.0, 224.0, 197.0, 194.0, 160.0, 151.0, 141.0, 90.0, 112.0, 102.0, 13.0, 125.0, 212.0, 60.0, 117.0, 131.0, 186.0, 258.0, 190.0, 165.0, 157.0, 146.0, 109.0, 30.0, 25.0, 37.0, 115.0, 183.0, 199.0, 210.0, 59.0, 281.0, 182.0, 156.0, 81.0, 60.0, 44.0, 43.0, 61.0, 106.0, 170.0, 198.0, 50.0, 36.0, 65.0, 133.0, 156.0, 95.0, 40.0, 64.0, 235.0, 176.0, 106.0, 44.0, 331.0, 175.0, 144.0, 75.0, 82.0, 61.0, 104.0, 144.0, 68.0, 108.0, 118.0, 360.0, 290.0, 95.0, 89.0, 85.0, 65.0, 43.0, 6.0, 30.0, 377.0, 359.0, 209.0, 116.0, 96.0, 70.0, 43.0, 291.0, 125.0, 118.0, 111.0, 65.0, 28.0, 259.0, 222.0, 182.0, 167.0, 142.0, 23.0, 15.0, 49.0, 373.0, 342.0, 141.0, 414.0, 317.0, 283.0, 167.0, 129.0, 43.0, 417.0, 332.0, 325.0, 132.0, 93.0, 416.0, 341.0, 241.0, 104.0, 33.0, 273.0, 144.0, 5.0, 338.0, 316.0, 276.0, 245.0, 230.0, 181.0, 166.0, 133.0, 93.0, 52.0, 35.0, 7.0, 435.0, 428.0, 418.0, 315.0, 296.0, 281.0, 242.0, 231.0, 122.0]
d = [1.0, 12.0, 9.0, 3.0, 2.0, 1.0, 18.0, 5.0, 7.0, 19.0, 16.0, 8.0, 19.0, 3.0, 3.0, 11.0, 2.0, 11.0, 3.0, 6.0, 14.0, 7.0, 4.0, 17.0, 13.0, 11.0, 3.0, 13.0, 7.0, 14.0, 11.0, 10.0, 19.0, 17.0, 6.0, 3.0, 15.0, 5.0, 17.0, 5.0, 9.0, 19.0, 5.0, 3.0, 6.0, 4.0, 17.0, 14.0, 6.0, 6.0, 19.0, 14.0, 11.0, 14.0, 16.0, 2.0, 10.0, 15.0, 15.0, 14.0, 17.0, 12.0, 11.0, 3.0, 1.0, 17.0, 19.0, 10.0, 18.0, 18.0, 12.0, 14.0, 17.0, 18.0, 9.0, 16.0, 4.0, 10.0, 4.0, 4.0, 2.0, 15.0, 6.0, 6.0, 14.0, 19.0, 4.0, 12.0, 10.0, 18.0, 16.0, 6.0, 3.0, 14.0, 9.0, 15.0, 3.0, 13.0, 20.0, 3.0, 5.0, 15.0, 8.0, 11.0, 15.0, 13.0, 11.0, 3.0, 20.0, 5.0, 17.0, 1.0, 17.0, 4.0, 3.0, 17.0, 16.0, 9.0, 11.0, 1.0, 1.0, 11.0, 8.0, 10.0, 15.0, 8.0, 16.0, 19.0, 12.0, 10.0, 7.0, 18.0, 4.0, 3.0, 1.0, 9.0, 4.0, 19.0, 7.0, 16.0, 6.0, 15.0, 13.0, 3.0, 8.0, 17.0, 5.0, 4.0, 6.0, 19.0, 20.0, 1.0, 13.0, 8.0, 10.0, 10.0, 14.0, 12.0, 13.0, 20.0, 3.0, 10.0, 2.0, 10.0, 15.0, 3.0, 7.0, 11.0, 4.0, 1.0, 7.0, 10.0, 2.0, 4.0, 14.0, 14.0, 2.0, 1.0, 8.0, 15.0, 9.0, 13.0, 12.0, 20.0, 19.0, 4.0, 4.0, 4.0, 12.0, 20.0, 19.0, 16.0, 14.0, 4.0, 10.0, 16.0, 9.0, 14.0, 12.0, 9.0, 19.0, 5.0, 7.0, 9.0, 1.0, 7.0, 11.0, 9.0, 6.0, 7.0, 5.0, 14.0, 17.0, 3.0, 9.0, 19.0, 19.0, 7.0, 10.0, 2.0, 1.0, 10.0, 15.0, 18.0, 13.0, 19.0, 14.0, 6.0, 14.0, 10.0, 20.0, 13.0, 2.0, 4.0, 14.0, 10.0, 15.0, 18.0, 2.0, 14.0, 5.0, 9.0, 6.0, 17.0, 13.0, 13.0, 7.0, 19.0, 19.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.826