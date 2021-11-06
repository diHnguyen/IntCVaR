edge = [1 4; 1 15; 1 22; 1 27; 1 29; 1 38; 1 41; 1 42; 1 43; 2 21; 2 27; 2 36; 2 39; 2 40; 3 7; 3 11; 3 15; 3 21; 3 32; 3 42; 4 11; 4 21; 4 39; 4 48; 4 50; 5 3; 5 10; 5 14; 5 17; 5 19; 6 4; 6 28; 7 17; 7 20; 7 23; 7 28; 7 45; 8 7; 8 34; 8 35; 8 36; 8 39; 8 43; 9 5; 9 10; 9 11; 9 15; 9 16; 9 17; 9 30; 9 33; 9 43; 9 46; 10 4; 10 13; 10 30; 10 36; 10 38; 10 44; 11 21; 11 25; 11 31; 11 42; 12 5; 12 13; 12 20; 12 36; 12 45; 13 6; 13 23; 13 32; 13 41; 13 45; 14 3; 14 11; 14 15; 14 22; 14 31; 15 16; 15 19; 16 18; 17 20; 17 22; 17 30; 17 47; 18 4; 18 11; 18 25; 18 32; 19 11; 19 12; 19 13; 20 10; 20 14; 20 18; 20 41; 20 46; 21 13; 21 15; 21 17; 21 40; 22 12; 22 15; 22 40; 22 43; 23 12; 23 16; 24 8; 24 14; 24 19; 24 28; 24 36; 24 42; 24 47; 25 12; 25 16; 25 23; 25 41; 26 21; 26 22; 26 44; 27 23; 27 24; 27 26; 27 31; 27 32; 27 45; 27 50; 28 26; 28 47; 28 48; 28 50; 29 13; 30 11; 30 22; 30 28; 30 36; 30 47; 31 15; 31 19; 31 38; 31 43; 31 50; 32 15; 32 19; 32 22; 32 23; 32 49; 33 13; 33 18; 33 35; 33 38; 33 40; 33 41; 33 45; 33 47; 34 17; 34 21; 34 24; 34 30; 35 5; 35 6; 35 27; 36 21; 36 33; 36 40; 36 44; 36 49; 37 8; 37 23; 37 30; 37 34; 37 38; 37 44; 37 48; 38 11; 38 19; 38 34; 38 50; 39 6; 39 29; 39 36; 40 25; 40 26; 40 35; 40 36; 40 48; 41 3; 41 6; 41 12; 41 18; 41 23; 41 39; 41 45; 42 20; 42 24; 42 36; 42 40; 42 41; 42 44; 43 12; 43 16; 43 17; 43 20; 43 41; 43 50; 44 15; 44 19; 44 29; 44 34; 45 23; 45 30; 45 39; 45 44; 46 9; 46 38; 47 26; 47 36; 47 37; 48 2; 48 17; 48 23; 48 25; 48 28; 48 40; 48 41; 48 42; 48 47; 48 49; 49 4; 49 10; 49 12; 49 47]
cL_orig = [31.0, 137.0, 205.0, 262.0, 280.0, 373.0, 396.0, 413.0, 415.0, 193.0, 249.0, 335.0, 369.0, 379.0, 41.0, 82.0, 124.0, 182.0, 288.0, 393.0, 70.0, 170.0, 346.0, 445.0, 465.0, 19.0, 48.0, 87.0, 118.0, 137.0, 16.0, 216.0, 102.0, 135.0, 164.0, 215.0, 384.0, 5.0, 255.0, 269.0, 285.0, 305.0, 352.0, 36.0, 8.0, 15.0, 58.0, 69.0, 76.0, 214.0, 245.0, 336.0, 368.0, 64.0, 33.0, 202.0, 256.0, 285.0, 335.0, 101.0, 144.0, 203.0, 314.0, 68.0, 13.0, 77.0, 237.0, 333.0, 67.0, 95.0, 191.0, 280.0, 322.0, 99.0, 33.0, 0.0, 81.0, 174.0, 12.0, 45.0, 22.0, 33.0, 44.0, 128.0, 304.0, 142.0, 72.0, 68.0, 140.0, 77.0, 65.0, 63.0, 97.0, 58.0, 21.0, 211.0, 246.0, 75.0, 59.0, 40.0, 195.0, 96.0, 68.0, 175.0, 211.0, 115.0, 68.0, 150.0, 97.0, 48.0, 44.0, 125.0, 181.0, 231.0, 131.0, 95.0, 20.0, 159.0, 50.0, 38.0, 177.0, 41.0, 31.0, 15.0, 37.0, 55.0, 185.0, 229.0, 16.0, 192.0, 197.0, 220.0, 158.0, 187.0, 79.0, 20.0, 60.0, 168.0, 149.0, 124.0, 70.0, 122.0, 192.0, 175.0, 133.0, 101.0, 94.0, 165.0, 202.0, 152.0, 20.0, 55.0, 67.0, 80.0, 120.0, 137.0, 160.0, 129.0, 99.0, 35.0, 295.0, 286.0, 73.0, 148.0, 26.0, 36.0, 78.0, 131.0, 285.0, 145.0, 74.0, 26.0, 15.0, 68.0, 115.0, 272.0, 193.0, 44.0, 119.0, 334.0, 95.0, 34.0, 149.0, 141.0, 49.0, 45.0, 75.0, 385.0, 346.0, 295.0, 232.0, 179.0, 25.0, 43.0, 219.0, 180.0, 58.0, 19.0, 3.0, 23.0, 308.0, 271.0, 247.0, 233.0, 20.0, 69.0, 286.0, 246.0, 144.0, 91.0, 220.0, 150.0, 61.0, 8.0, 375.0, 77.0, 208.0, 110.0, 95.0, 462.0, 314.0, 251.0, 226.0, 198.0, 84.0, 70.0, 58.0, 8.0, 9.0, 445.0, 386.0, 371.0, 20.0]
cU_orig = [31.0, 137.0, 205.0, 262.0, 280.0, 373.0, 396.0, 413.0, 415.0, 193.0, 249.0, 335.0, 369.0, 379.0, 41.0, 82.0, 124.0, 182.0, 288.0, 393.0, 70.0, 170.0, 346.0, 445.0, 465.0, 19.0, 60.0, 99.0, 118.0, 137.0, 32.0, 216.0, 102.0, 135.0, 164.0, 215.0, 384.0, 5.0, 255.0, 269.0, 285.0, 305.0, 352.0, 36.0, 8.0, 15.0, 58.0, 69.0, 76.0, 214.0, 245.0, 350.0, 368.0, 64.0, 33.0, 202.0, 256.0, 285.0, 351.0, 101.0, 144.0, 203.0, 314.0, 68.0, 13.0, 89.0, 237.0, 333.0, 67.0, 109.0, 191.0, 280.0, 322.0, 111.0, 33.0, 18.0, 81.0, 174.0, 12.0, 45.0, 22.0, 33.0, 64.0, 128.0, 304.0, 142.0, 72.0, 68.0, 140.0, 77.0, 65.0, 63.0, 97.0, 58.0, 21.0, 211.0, 266.0, 75.0, 59.0, 40.0, 195.0, 96.0, 68.0, 175.0, 211.0, 115.0, 68.0, 162.0, 97.0, 48.0, 44.0, 125.0, 181.0, 231.0, 131.0, 95.0, 20.0, 159.0, 50.0, 38.0, 177.0, 41.0, 31.0, 15.0, 37.0, 55.0, 185.0, 229.0, 16.0, 192.0, 197.0, 220.0, 158.0, 187.0, 79.0, 20.0, 60.0, 168.0, 161.0, 124.0, 70.0, 122.0, 192.0, 175.0, 133.0, 101.0, 94.0, 165.0, 202.0, 152.0, 20.0, 55.0, 67.0, 80.0, 120.0, 137.0, 174.0, 129.0, 99.0, 35.0, 295.0, 286.0, 87.0, 148.0, 44.0, 36.0, 78.0, 131.0, 285.0, 145.0, 74.0, 26.0, 15.0, 78.0, 115.0, 272.0, 193.0, 44.0, 119.0, 334.0, 95.0, 34.0, 149.0, 141.0, 49.0, 45.0, 75.0, 385.0, 346.0, 295.0, 232.0, 179.0, 25.0, 43.0, 219.0, 180.0, 58.0, 19.0, 23.0, 23.0, 308.0, 271.0, 265.0, 233.0, 20.0, 69.0, 296.0, 246.0, 156.0, 103.0, 220.0, 150.0, 61.0, 8.0, 375.0, 77.0, 208.0, 110.0, 95.0, 462.0, 314.0, 251.0, 226.0, 212.0, 84.0, 70.0, 58.0, 8.0, 9.0, 445.0, 386.0, 371.0, 20.0]
d = [4.0, 6.0, 18.0, 14.0, 6.0, 10.0, 7.0, 15.0, 15.0, 12.0, 4.0, 1.0, 1.0, 9.0, 11.0, 2.0, 15.0, 15.0, 11.0, 4.0, 17.0, 7.0, 16.0, 9.0, 1.0, 17.0, 13.0, 3.0, 6.0, 9.0, 9.0, 14.0, 19.0, 20.0, 16.0, 18.0, 1.0, 10.0, 16.0, 1.0, 4.0, 18.0, 10.0, 15.0, 1.0, 18.0, 18.0, 13.0, 2.0, 9.0, 3.0, 6.0, 19.0, 10.0, 3.0, 13.0, 17.0, 6.0, 7.0, 15.0, 18.0, 19.0, 2.0, 17.0, 4.0, 6.0, 3.0, 6.0, 7.0, 10.0, 15.0, 10.0, 13.0, 14.0, 19.0, 20.0, 15.0, 11.0, 7.0, 15.0, 13.0, 17.0, 4.0, 11.0, 9.0, 9.0, 6.0, 4.0, 20.0, 20.0, 4.0, 3.0, 10.0, 19.0, 8.0, 19.0, 7.0, 9.0, 6.0, 1.0, 19.0, 4.0, 18.0, 20.0, 17.0, 7.0, 10.0, 15.0, 20.0, 15.0, 7.0, 20.0, 10.0, 12.0, 20.0, 15.0, 2.0, 20.0, 12.0, 10.0, 4.0, 9.0, 15.0, 18.0, 19.0, 7.0, 20.0, 20.0, 15.0, 16.0, 15.0, 16.0, 9.0, 4.0, 3.0, 18.0, 9.0, 7.0, 4.0, 9.0, 5.0, 2.0, 18.0, 20.0, 1.0, 14.0, 12.0, 17.0, 13.0, 6.0, 3.0, 17.0, 7.0, 8.0, 16.0, 1.0, 5.0, 12.0, 4.0, 1.0, 12.0, 9.0, 16.0, 4.0, 10.0, 17.0, 20.0, 17.0, 6.0, 15.0, 13.0, 8.0, 11.0, 18.0, 8.0, 10.0, 20.0, 19.0, 15.0, 12.0, 15.0, 16.0, 3.0, 8.0, 6.0, 9.0, 18.0, 10.0, 7.0, 8.0, 6.0, 8.0, 18.0, 10.0, 10.0, 8.0, 13.0, 16.0, 20.0, 17.0, 9.0, 13.0, 18.0, 3.0, 15.0, 2.0, 2.0, 2.0, 15.0, 12.0, 7.0, 10.0, 6.0, 9.0, 13.0, 6.0, 20.0, 3.0, 19.0, 19.0, 3.0, 1.0, 13.0, 20.0, 19.0, 12.0, 5.0, 8.0, 8.0, 9.0, 2.0, 16.0, 12.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.408
