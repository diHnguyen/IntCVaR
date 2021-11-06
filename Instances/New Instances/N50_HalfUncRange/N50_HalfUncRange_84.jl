edge = [1 6; 1 13; 1 14; 1 27; 1 33; 1 42; 2 7; 2 17; 2 35; 2 37; 2 41; 2 43; 3 5; 3 7; 3 14; 3 49; 4 2; 4 26; 4 31; 4 35; 4 40; 5 35; 5 37; 5 39; 6 19; 6 25; 6 37; 6 41; 6 45; 7 4; 7 5; 7 15; 7 23; 7 39; 7 41; 7 42; 7 48; 8 17; 8 18; 8 29; 9 5; 9 7; 9 10; 9 13; 9 14; 9 20; 9 21; 9 25; 9 31; 9 32; 9 39; 9 43; 9 50; 10 6; 10 27; 10 28; 10 31; 10 40; 10 44; 10 49; 11 14; 11 25; 11 42; 11 46; 12 20; 12 23; 12 33; 13 6; 13 12; 13 14; 13 17; 13 31; 13 38; 13 40; 13 42; 14 5; 14 16; 14 23; 14 24; 14 27; 14 35; 15 2; 15 9; 15 21; 15 32; 15 36; 15 43; 16 10; 16 20; 16 37; 17 18; 17 22; 17 27; 17 43; 17 49; 18 14; 18 17; 18 21; 18 27; 18 42; 18 50; 19 5; 19 13; 19 20; 19 29; 20 6; 20 15; 20 21; 20 30; 20 46; 21 3; 21 4; 21 5; 21 29; 21 38; 21 46; 22 24; 22 34; 22 48; 23 14; 23 40; 23 41; 23 46; 24 13; 24 18; 24 45; 24 47; 25 8; 25 10; 25 13; 25 15; 25 21; 25 23; 25 30; 26 14; 26 15; 26 25; 26 31; 26 38; 26 50; 27 30; 27 36; 27 41; 27 49; 28 4; 28 8; 28 12; 28 18; 28 30; 28 34; 28 36; 28 42; 28 45; 29 5; 29 12; 29 18; 29 20; 29 21; 29 50; 30 8; 30 10; 30 21; 30 23; 30 25; 30 29; 30 33; 30 37; 31 8; 31 13; 31 25; 31 44; 31 50; 32 2; 32 18; 32 24; 32 37; 32 43; 32 45; 32 46; 33 6; 33 9; 33 10; 33 31; 33 41; 33 46; 33 47; 34 22; 34 26; 34 40; 35 9; 35 19; 35 24; 35 26; 36 8; 36 13; 36 25; 36 26; 36 28; 36 44; 36 47; 37 11; 37 23; 37 24; 37 32; 37 36; 38 7; 38 10; 38 19; 38 23; 38 34; 38 39; 38 43; 39 9; 39 12; 39 16; 39 34; 39 38; 39 46; 39 48; 39 50; 40 3; 40 10; 40 15; 40 28; 40 32; 40 34; 40 37; 40 42; 40 46; 41 3; 41 4; 41 15; 41 18; 41 21; 42 9; 42 30; 43 28; 43 29; 43 32; 43 40; 44 2; 44 4; 44 7; 44 11; 44 16; 44 28; 44 32; 45 10; 45 23; 45 35; 45 41; 45 42; 46 14; 46 18; 46 30; 46 34; 47 5; 47 21; 47 33; 47 44; 47 48; 48 3; 48 21; 48 24; 48 30; 49 9; 49 15; 49 21; 49 23; 49 25; 49 26; 49 28; 49 37; 49 39]
cL_orig = [45.0, 114.0, 129.0, 261.0, 319.0, 413.0, 49.0, 150.0, 335.0, 354.0, 384.0, 410.0, 17.0, 42.0, 107.0, 461.0, 17.0, 213.5, 268.0, 305.0, 358.0, 297.0, 318.0, 338.0, 132.0, 183.0, 309.0, 353.0, 389.0, 28.0, 13.5, 73.0, 160.0, 310.0, 337.0, 353.0, 406.0, 93.0, 104.0, 210.0, 43.0, 18.0, 12.0, 44.0, 49.5, 105.0, 118.0, 160.0, 216.0, 235.0, 297.0, 341.0, 415.0, 39.0, 175.0, 180.5, 212.0, 305.0, 340.0, 383.0, 26.0, 144.0, 307.0, 353.0, 75.0, 108.0, 211.0, 75.0, 8.0, 9.0, 43.0, 171.5, 255.0, 261.0, 295.0, 92.0, 17.0, 88.0, 103.0, 128.5, 215.0, 133.0, 62.0, 62.0, 166.0, 211.0, 282.0, 60.0, 39.0, 209.0, 15.0, 53.0, 97.0, 259.0, 323.0, 41.0, 10.5, 26.0, 90.0, 240.0, 324.0, 144.0, 56.0, 7.0, 96.5, 140.0, 55.0, 5.0, 98.0, 261.0, 181.0, 171.0, 158.0, 84.0, 166.0, 246.0, 19.0, 123.0, 261.5, 92.0, 167.0, 184.0, 229.0, 113.0, 60.0, 202.0, 226.0, 170.0, 145.0, 125.0, 96.0, 36.0, 14.5, 49.5, 119.0, 108.0, 5.0, 54.0, 119.0, 236.5, 27.0, 85.0, 143.0, 220.0, 235.0, 196.0, 157.0, 103.0, 21.5, 59.0, 77.0, 136.0, 165.0, 240.0, 167.5, 110.0, 92.0, 82.0, 202.0, 219.5, 204.0, 89.0, 75.0, 46.0, 5.0, 26.0, 72.0, 226.0, 185.0, 58.0, 135.0, 194.0, 305.0, 137.0, 72.5, 52.0, 105.0, 135.0, 136.0, 266.0, 241.0, 228.0, 17.0, 83.0, 126.0, 136.0, 123.0, 73.0, 63.0, 263.0, 162.0, 113.0, 85.0, 275.0, 230.0, 110.0, 101.0, 85.0, 85.0, 105.0, 260.0, 130.0, 130.5, 49.0, 15.0, 306.0, 280.0, 191.0, 151.0, 37.0, 8.0, 45.0, 297.0, 273.0, 226.0, 50.0, 8.0, 70.0, 86.0, 109.0, 372.0, 302.0, 246.0, 122.0, 84.0, 57.0, 31.0, 19.0, 60.0, 379.0, 369.0, 262.0, 235.0, 205.0, 328.0, 117.0, 146.0, 139.0, 107.0, 35.0, 412.0, 400.0, 362.0, 331.0, 279.0, 151.0, 116.0, 354.0, 216.0, 102.0, 40.0, 30.0, 318.0, 280.0, 158.0, 111.5, 425.0, 259.0, 140.0, 29.0, 11.0, 454.0, 268.0, 244.0, 185.0, 403.0, 338.0, 279.0, 257.0, 237.5, 223.5, 202.5, 119.0, 99.0]
cU_orig = [45.0, 122.0, 129.0, 261.0, 319.0, 413.0, 55.0, 150.0, 335.0, 354.0, 392.0, 410.0, 17.0, 42.0, 107.0, 461.0, 17.0, 222.5, 268.0, 305.0, 358.0, 297.0, 318.0, 338.0, 132.0, 193.0, 319.0, 353.0, 389.0, 28.0, 18.5, 79.0, 160.0, 320.0, 337.0, 353.0, 406.0, 93.0, 104.0, 220.0, 43.0, 18.0, 12.0, 44.0, 56.5, 105.0, 118.0, 168.0, 216.0, 235.0, 297.0, 341.0, 415.0, 39.0, 175.0, 185.5, 212.0, 305.0, 340.0, 389.0, 26.0, 144.0, 307.0, 353.0, 75.0, 108.0, 211.0, 75.0, 8.0, 9.0, 43.0, 178.5, 255.0, 271.0, 295.0, 92.0, 17.0, 88.0, 103.0, 135.5, 215.0, 133.0, 62.0, 62.0, 166.0, 211.0, 282.0, 60.0, 39.0, 209.0, 15.0, 53.0, 97.0, 259.0, 323.0, 41.0, 17.5, 26.0, 98.0, 240.0, 324.0, 144.0, 56.0, 7.0, 105.5, 140.0, 55.0, 15.0, 98.0, 261.0, 181.0, 171.0, 158.0, 84.0, 166.0, 256.0, 19.0, 123.0, 266.5, 98.0, 167.0, 184.0, 229.0, 113.0, 60.0, 210.0, 226.0, 170.0, 145.0, 125.0, 96.0, 36.0, 21.5, 58.5, 119.0, 108.0, 5.0, 54.0, 119.0, 245.5, 27.0, 85.0, 143.0, 220.0, 235.0, 196.0, 157.0, 103.0, 28.5, 59.0, 77.0, 136.0, 165.0, 240.0, 176.5, 116.0, 92.0, 82.0, 212.0, 228.5, 204.0, 89.0, 75.0, 56.0, 5.0, 26.0, 72.0, 234.0, 185.0, 68.0, 135.0, 194.0, 305.0, 143.0, 81.5, 52.0, 105.0, 135.0, 136.0, 266.0, 241.0, 238.0, 17.0, 83.0, 126.0, 136.0, 123.0, 79.0, 63.0, 263.0, 162.0, 113.0, 85.0, 275.0, 230.0, 110.0, 101.0, 85.0, 85.0, 105.0, 260.0, 140.0, 137.5, 49.0, 15.0, 306.0, 280.0, 191.0, 151.0, 37.0, 8.0, 45.0, 305.0, 273.0, 226.0, 50.0, 8.0, 70.0, 86.0, 109.0, 378.0, 302.0, 246.0, 122.0, 84.0, 57.0, 31.0, 19.0, 60.0, 379.0, 379.0, 262.0, 235.0, 205.0, 328.0, 117.0, 154.0, 139.0, 107.0, 35.0, 420.0, 400.0, 372.0, 331.0, 287.0, 161.0, 116.0, 354.0, 216.0, 102.0, 40.0, 30.0, 318.0, 280.0, 168.0, 118.5, 425.0, 259.0, 140.0, 29.0, 11.0, 454.0, 268.0, 244.0, 185.0, 403.0, 338.0, 279.0, 257.0, 244.5, 228.5, 209.5, 119.0, 99.0]
d = [4.0, 1.0, 5.0, 5.0, 12.0, 16.0, 14.0, 5.0, 20.0, 14.0, 5.0, 15.0, 19.0, 2.0, 4.0, 20.0, 14.0, 11.0, 11.0, 18.0, 7.0, 14.0, 7.0, 15.0, 3.0, 15.0, 17.0, 10.0, 7.0, 2.0, 5.0, 3.0, 16.0, 19.0, 18.0, 9.0, 16.0, 18.0, 13.0, 19.0, 16.0, 19.0, 2.0, 16.0, 14.0, 11.0, 20.0, 5.0, 16.0, 1.0, 9.0, 12.0, 5.0, 7.0, 5.0, 15.0, 8.0, 19.0, 5.0, 11.0, 13.0, 14.0, 5.0, 8.0, 12.0, 12.0, 11.0, 9.0, 20.0, 20.0, 14.0, 11.0, 2.0, 10.0, 2.0, 5.0, 2.0, 6.0, 17.0, 6.0, 4.0, 13.0, 20.0, 12.0, 15.0, 19.0, 14.0, 4.0, 4.0, 7.0, 19.0, 8.0, 8.0, 16.0, 14.0, 19.0, 17.0, 9.0, 15.0, 3.0, 12.0, 19.0, 17.0, 14.0, 9.0, 2.0, 16.0, 14.0, 18.0, 1.0, 4.0, 5.0, 20.0, 17.0, 17.0, 14.0, 6.0, 20.0, 17.0, 19.0, 12.0, 17.0, 1.0, 17.0, 7.0, 9.0, 11.0, 20.0, 16.0, 14.0, 14.0, 1.0, 16.0, 10.0, 18.0, 7.0, 15.0, 12.0, 17.0, 12.0, 9.0, 11.0, 19.0, 7.0, 6.0, 17.0, 11.0, 19.0, 3.0, 15.0, 8.0, 8.0, 4.0, 19.0, 10.0, 6.0, 15.0, 13.0, 7.0, 1.0, 6.0, 11.0, 12.0, 20.0, 12.0, 8.0, 15.0, 20.0, 5.0, 2.0, 4.0, 5.0, 10.0, 1.0, 2.0, 9.0, 8.0, 19.0, 18.0, 13.0, 8.0, 1.0, 3.0, 19.0, 20.0, 15.0, 19.0, 16.0, 10.0, 3.0, 10.0, 1.0, 7.0, 12.0, 2.0, 20.0, 16.0, 9.0, 11.0, 19.0, 19.0, 17.0, 3.0, 20.0, 5.0, 7.0, 10.0, 13.0, 12.0, 16.0, 1.0, 20.0, 4.0, 5.0, 9.0, 5.0, 11.0, 5.0, 10.0, 19.0, 17.0, 10.0, 7.0, 9.0, 12.0, 4.0, 14.0, 12.0, 18.0, 4.0, 13.0, 7.0, 17.0, 9.0, 3.0, 15.0, 6.0, 14.0, 13.0, 3.0, 3.0, 18.0, 10.0, 16.0, 13.0, 13.0, 8.0, 7.0, 6.0, 15.0, 12.0, 20.0, 5.0, 8.0, 1.0, 8.0, 4.0, 16.0, 6.0, 12.0, 17.0, 11.0, 5.0, 20.0, 12.0, 14.0, 17.0, 13.0, 9.0, 19.0, 4.0, 6.0, 1.0, 1.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
β = 0.04
