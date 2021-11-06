edge = [1 6; 1 10; 1 22; 1 27; 1 34; 1 39; 1 43; 1 45; 1 50; 2 4; 2 6; 2 18; 2 27; 2 29; 2 30; 2 40; 2 45; 3 5; 3 6; 3 11; 3 13; 3 15; 3 21; 3 25; 3 34; 4 25; 4 30; 4 35; 4 39; 4 40; 4 45; 4 46; 5 10; 5 12; 5 20; 5 21; 5 25; 5 31; 5 44; 5 46; 6 10; 6 12; 6 13; 6 14; 6 21; 6 44; 6 45; 7 12; 7 15; 7 27; 7 33; 7 41; 7 42; 7 47; 8 10; 8 17; 8 18; 8 19; 8 20; 8 23; 8 40; 8 46; 8 50; 9 37; 10 7; 10 20; 10 22; 10 29; 10 30; 11 32; 11 40; 11 41; 11 42; 11 44; 12 7; 12 40; 12 41; 13 10; 13 22; 13 24; 13 40; 14 21; 14 33; 14 39; 15 2; 15 19; 15 24; 15 30; 16 18; 16 26; 16 31; 16 34; 16 37; 17 22; 17 30; 17 33; 17 47; 18 14; 18 35; 18 41; 18 42; 19 10; 19 23; 19 30; 19 34; 19 42; 20 7; 20 15; 20 37; 20 43; 20 45; 21 6; 21 18; 21 28; 21 31; 22 13; 22 21; 22 31; 22 35; 22 37; 23 3; 23 7; 23 27; 23 43; 24 14; 25 23; 25 34; 25 42; 25 46; 26 11; 26 25; 26 29; 26 41; 27 21; 27 22; 27 33; 27 36; 27 43; 28 4; 28 11; 28 23; 28 26; 28 31; 29 4; 29 6; 29 8; 29 37; 29 45; 30 2; 30 9; 30 16; 30 28; 30 34; 30 42; 30 43; 30 47; 30 49; 30 50; 31 11; 32 5; 32 18; 32 39; 33 17; 33 19; 33 20; 33 25; 33 42; 34 3; 34 10; 35 8; 36 5; 36 12; 36 13; 36 22; 36 23; 36 39; 36 41; 36 46; 36 47; 37 14; 37 28; 37 31; 37 41; 38 3; 38 19; 38 20; 38 28; 38 34; 38 47; 39 7; 39 15; 39 19; 39 27; 39 34; 39 35; 39 48; 40 49; 40 50; 41 44; 42 13; 42 16; 42 39; 43 15; 43 21; 43 28; 43 35; 43 41; 43 44; 44 5; 44 10; 44 34; 44 38; 44 39; 44 40; 45 16; 45 19; 45 20; 45 21; 45 29; 45 34; 45 41; 45 44; 46 16; 46 20; 46 21; 46 37; 47 3; 47 5; 47 10; 47 21; 48 3; 48 8; 48 23; 48 41; 49 8; 49 31; 49 39]
cL_orig = [49.0, 88.0, 208.0, 262.0, 333.0, 377.0, 420.0, 443.0, 485.0, 22.0, 44.0, 158.0, 247.0, 268.0, 284.0, 377.0, 427.0, 19.0, 29.0, 79.0, 99.0, 125.0, 180.0, 217.0, 311.0, 205.0, 262.0, 313.0, 347.0, 362.0, 412.0, 424.0, 48.0, 73.0, 152.0, 157.0, 201.0, 261.0, 385.0, 415.0, 41.0, 57.0, 70.0, 76.0, 147.0, 373.0, 390.0, 52.0, 78.0, 199.0, 260.0, 337.0, 350.0, 402.0, 19.0, 87.0, 104.0, 114.0, 117.0, 148.0, 319.0, 382.0, 425.0, 284.0, 26.0, 103.0, 117.0, 185.0, 198.0, 207.0, 285.0, 296.0, 312.0, 334.0, 51.0, 285.0, 293.0, 35.0, 89.0, 109.0, 260.0, 70.0, 185.0, 255.0, 128.0, 41.0, 94.0, 150.0, 7.0, 104.0, 150.0, 182.0, 212.0, 54.0, 132.0, 160.0, 295.0, 32.0, 173.0, 220.0, 244.0, 88.0, 41.0, 110.0, 154.0, 232.0, 127.0, 55.0, 159.0, 226.0, 248.0, 151.0, 27.0, 70.0, 100.0, 95.0, 11.0, 93.0, 134.0, 145.0, 199.0, 164.0, 35.0, 196.0, 100.0, 13.0, 79.0, 174.0, 214.0, 148.0, 12.0, 32.0, 150.0, 60.0, 54.0, 51.0, 95.0, 155.0, 228.0, 171.0, 49.0, 20.0, 26.0, 247.0, 235.0, 215.0, 78.0, 165.0, 280.0, 212.0, 145.0, 15.0, 35.0, 120.0, 134.0, 174.0, 193.0, 198.0, 198.0, 266.0, 138.0, 72.0, 161.0, 126.0, 128.0, 79.0, 91.0, 307.0, 237.0, 260.0, 310.0, 237.0, 222.0, 136.0, 117.0, 27.0, 54.0, 89.0, 107.0, 234.0, 92.0, 60.0, 39.0, 345.0, 180.0, 177.0, 105.0, 38.0, 83.0, 324.0, 235.0, 203.0, 115.0, 45.0, 36.0, 95.0, 86.0, 104.0, 32.0, 290.0, 263.0, 32.0, 275.0, 221.0, 154.0, 82.0, 20.0, 5.0, 394.0, 339.0, 105.0, 63.0, 53.0, 43.0, 292.0, 256.0, 247.0, 244.0, 164.0, 106.0, 44.0, 12.0, 292.0, 258.0, 253.0, 92.0, 442.0, 421.0, 368.0, 257.0, 441.0, 397.0, 246.0, 72.0, 408.0, 185.0, 98.0]
cU_orig = [49.0, 88.0, 208.0, 262.0, 333.0, 377.0, 420.0, 443.0, 503.0, 22.0, 44.0, 158.0, 247.0, 268.0, 284.0, 377.0, 443.0, 19.0, 29.0, 79.0, 99.0, 125.0, 180.0, 217.0, 311.0, 205.0, 262.0, 313.0, 347.0, 362.0, 412.0, 424.0, 48.0, 73.0, 152.0, 157.0, 201.0, 261.0, 385.0, 415.0, 41.0, 57.0, 70.0, 76.0, 147.0, 393.0, 390.0, 52.0, 78.0, 211.0, 260.0, 337.0, 350.0, 402.0, 19.0, 87.0, 104.0, 114.0, 117.0, 148.0, 319.0, 382.0, 425.0, 284.0, 26.0, 103.0, 117.0, 201.0, 198.0, 207.0, 285.0, 296.0, 312.0, 334.0, 51.0, 285.0, 293.0, 35.0, 89.0, 109.0, 276.0, 70.0, 185.0, 255.0, 128.0, 41.0, 94.0, 150.0, 25.0, 104.0, 150.0, 182.0, 212.0, 54.0, 132.0, 160.0, 295.0, 46.0, 173.0, 234.0, 244.0, 88.0, 41.0, 110.0, 154.0, 232.0, 127.0, 55.0, 173.0, 226.0, 248.0, 151.0, 27.0, 70.0, 100.0, 95.0, 11.0, 93.0, 134.0, 145.0, 199.0, 164.0, 35.0, 196.0, 100.0, 23.0, 95.0, 174.0, 214.0, 148.0, 12.0, 32.0, 150.0, 60.0, 54.0, 63.0, 95.0, 155.0, 242.0, 171.0, 49.0, 20.0, 26.0, 247.0, 235.0, 215.0, 78.0, 165.0, 280.0, 212.0, 145.0, 15.0, 35.0, 120.0, 134.0, 174.0, 193.0, 198.0, 198.0, 266.0, 138.0, 72.0, 161.0, 146.0, 128.0, 79.0, 91.0, 307.0, 237.0, 278.0, 310.0, 237.0, 232.0, 136.0, 135.0, 27.0, 54.0, 107.0, 107.0, 234.0, 92.0, 60.0, 39.0, 345.0, 190.0, 177.0, 105.0, 38.0, 95.0, 324.0, 235.0, 203.0, 115.0, 45.0, 36.0, 95.0, 86.0, 104.0, 32.0, 290.0, 263.0, 32.0, 275.0, 221.0, 154.0, 82.0, 20.0, 17.0, 394.0, 339.0, 105.0, 63.0, 53.0, 43.0, 292.0, 256.0, 247.0, 244.0, 164.0, 106.0, 44.0, 12.0, 306.0, 258.0, 253.0, 92.0, 442.0, 421.0, 368.0, 257.0, 455.0, 397.0, 246.0, 72.0, 408.0, 185.0, 98.0]
d = [14.0, 3.0, 5.0, 16.0, 6.0, 11.0, 13.0, 14.0, 7.0, 14.0, 11.0, 19.0, 3.0, 10.0, 5.0, 17.0, 4.0, 16.0, 10.0, 17.0, 11.0, 2.0, 11.0, 5.0, 4.0, 3.0, 18.0, 19.0, 2.0, 13.0, 13.0, 16.0, 19.0, 2.0, 6.0, 8.0, 10.0, 4.0, 19.0, 16.0, 17.0, 10.0, 14.0, 20.0, 10.0, 2.0, 2.0, 10.0, 13.0, 8.0, 11.0, 6.0, 3.0, 20.0, 8.0, 20.0, 6.0, 5.0, 2.0, 2.0, 16.0, 7.0, 3.0, 5.0, 2.0, 7.0, 6.0, 8.0, 9.0, 3.0, 5.0, 12.0, 9.0, 8.0, 3.0, 19.0, 6.0, 17.0, 3.0, 9.0, 18.0, 13.0, 16.0, 11.0, 16.0, 11.0, 20.0, 3.0, 17.0, 13.0, 6.0, 5.0, 15.0, 2.0, 14.0, 14.0, 14.0, 10.0, 19.0, 5.0, 15.0, 12.0, 19.0, 10.0, 7.0, 18.0, 5.0, 6.0, 11.0, 18.0, 15.0, 14.0, 4.0, 3.0, 6.0, 13.0, 6.0, 14.0, 20.0, 19.0, 5.0, 5.0, 10.0, 9.0, 7.0, 15.0, 14.0, 1.0, 8.0, 5.0, 6.0, 9.0, 4.0, 10.0, 17.0, 17.0, 19.0, 13.0, 2.0, 3.0, 6.0, 19.0, 1.0, 18.0, 10.0, 3.0, 12.0, 20.0, 16.0, 18.0, 4.0, 10.0, 18.0, 4.0, 6.0, 12.0, 12.0, 17.0, 10.0, 5.0, 3.0, 1.0, 12.0, 4.0, 2.0, 15.0, 10.0, 4.0, 9.0, 5.0, 4.0, 19.0, 12.0, 9.0, 9.0, 16.0, 14.0, 3.0, 11.0, 8.0, 12.0, 15.0, 1.0, 19.0, 2.0, 1.0, 7.0, 8.0, 14.0, 13.0, 2.0, 7.0, 11.0, 18.0, 2.0, 4.0, 12.0, 15.0, 14.0, 14.0, 13.0, 11.0, 20.0, 12.0, 7.0, 11.0, 3.0, 19.0, 9.0, 10.0, 14.0, 14.0, 17.0, 8.0, 14.0, 18.0, 11.0, 4.0, 13.0, 8.0, 2.0, 4.0, 1.0, 19.0, 7.0, 8.0, 1.0, 14.0, 13.0, 5.0, 6.0, 13.0, 8.0, 7.0, 1.0, 2.0, 16.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.493