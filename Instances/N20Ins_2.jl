edge = [1 6; 1 13; 1 14; 1 18; 2 3; 2 6; 2 7; 2 8; 2 10; 2 11; 2 18; 2 20; 3 2; 3 4; 3 12; 4 3; 4 8; 4 16; 4 20; 5 6; 5 9; 5 10; 5 12; 5 13; 5 15; 5 16; 5 17; 5 19; 5 20; 6 9; 6 11; 6 13; 6 19; 7 2; 7 3; 7 5; 7 6; 7 10; 7 13; 7 15; 7 16; 7 17; 7 18; 7 19; 8 2; 8 6; 8 9; 8 10; 8 12; 8 17; 8 20; 9 6; 9 8; 9 11; 9 13; 9 14; 9 15; 9 17; 9 20; 10 3; 10 5; 10 6; 10 12; 10 13; 10 14; 10 16; 10 17; 10 20; 11 7; 11 8; 11 14; 11 18; 11 19; 12 2; 12 4; 12 14; 12 16; 12 20; 13 2; 13 3; 13 4; 13 5; 13 7; 13 11; 13 16; 13 20; 14 3; 14 8; 14 9; 14 12; 14 18; 14 19; 14 20; 15 3; 15 5; 15 6; 15 10; 15 17; 15 18; 16 6; 16 7; 16 9; 16 10; 16 19; 16 20; 17 2; 17 4; 17 5; 17 8; 17 9; 17 10; 17 11; 17 13; 17 14; 17 18; 18 4; 18 7; 18 15; 19 3; 19 4; 19 9; 19 14; 19 16; 19 18]
cL_orig = [8.0, 55.0, 116.0, 31.0, 41.0, 30.0, 178.0, 70.0, 254.0, 356.0, 165.0, 54.0, 43.0, 29.0, 440.0, 8.0, 68.0, 365.0, 772.0, 10.0, 181.0, 0.0, 8.0, 218.0, 335.0, 335.0, 9.0, 539.0, 523.0, 7.0, 205.0, 118.0, 337.0, 225.0, 3.0, 0.0, 10.0, 18.0, 48.0, 346.0, 94.0, 498.0, 124.0, 156.0, 214.0, 66.0, 6.0, 5.0, 54.0, 26.0, 222.0, 6.0, 40.0, 48.0, 54.0, 47.0, 59.0, 20.0, 341.0, 319.0, 52.0, 119.0, 11.0, 21.0, 77.0, 253.0, 312.0, 150.0, 89.0, 32.0, 89.0, 38.0, 11.0, 338.0, 247.0, 91.0, 58.0, 223.0, 516.0, 336.0, 113.0, 123.0, 9.0, 67.0, 22.0, 299.0, 148.0, 169.0, 148.0, 81.0, 68.0, 88.0, 22.0, 25.0, 344.0, 88.0, 89.0, 51.0, 10.0, 103.0, 437.0, 330.0, 208.0, 54.0, 68.0, 745.0, 212.0, 577.0, 189.0, 58.0, 346.0, 70.0, 19.0, 71.0, 34.0, 589.0, 40.0, 87.0, 530.0, 119.0, 22.0, 148.0, 21.0, 19.0]
cU_orig = [152.0, 389.0, 116.0, 1069.0, 41.0, 42.0, 178.0, 70.0, 254.0, 356.0, 165.0, 184.0, 57.0, 29.0, 458.0, 8.0, 68.0, 365.0, 772.0, 10.0, 181.0, 192.0, 348.0, 218.0, 335.0, 335.0, 9.0, 539.0, 523.0, 33.0, 205.0, 118.0, 337.0, 225.0, 59.0, 4.0, 10.0, 18.0, 474.0, 346.0, 466.0, 498.0, 306.0, 156.0, 214.0, 66.0, 18.0, 5.0, 54.0, 444.0, 440.0, 6.0, 40.0, 48.0, 54.0, 47.0, 59.0, 20.0, 341.0, 319.0, 146.0, 119.0, 11.0, 97.0, 77.0, 253.0, 312.0, 150.0, 227.0, 32.0, 89.0, 174.0, 11.0, 338.0, 509.0, 91.0, 58.0, 223.0, 516.0, 336.0, 113.0, 123.0, 9.0, 67.0, 22.0, 299.0, 848.0, 169.0, 250.0, 81.0, 68.0, 88.0, 22.0, 25.0, 344.0, 88.0, 267.0, 139.0, 10.0, 103.0, 437.0, 330.0, 208.0, 54.0, 270.0, 745.0, 1058.0, 577.0, 205.0, 234.0, 346.0, 70.0, 19.0, 71.0, 34.0, 589.0, 40.0, 87.0, 530.0, 119.0, 114.0, 148.0, 27.0, 21.0]
d = [30.0, 13.0, 33.0, 42.0, 47.0, 47.0, 38.0, 42.0, 40.0, 16.0, 39.0, 40.0, 35.0, 20.0, 33.0, 1.0, 43.0, 3.0, 41.0, 6.0, 45.0, 48.0, 37.0, 17.0, 3.0, 39.0, 31.0, 6.0, 11.0, 40.0, 27.0, 46.0, 14.0, 7.0, 10.0, 39.0, 37.0, 34.0, 16.0, 35.0, 22.0, 8.0, 9.0, 28.0, 35.0, 32.0, 24.0, 1.0, 1.0, 49.0, 18.0, 23.0, 27.0, 6.0, 9.0, 38.0, 36.0, 37.0, 32.0, 25.0, 19.0, 23.0, 49.0, 25.0, 34.0, 24.0, 5.0, 46.0, 40.0, 22.0, 16.0, 30.0, 1.0, 35.0, 10.0, 43.0, 19.0, 39.0, 22.0, 14.0, 14.0, 17.0, 24.0, 34.0, 50.0, 34.0, 16.0, 3.0, 3.0, 37.0, 24.0, 8.0, 49.0, 39.0, 27.0, 42.0, 6.0, 38.0, 28.0, 18.0, 46.0, 35.0, 24.0, 36.0, 43.0, 37.0, 18.0, 40.0, 31.0, 6.0, 49.0, 25.0, 11.0, 48.0, 23.0, 42.0, 6.0, 26.0, 34.0, 32.0, 27.0, 32.0, 17.0, 33.0]
Len = length(d)

yy = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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
delta2 = 2
last_node = maximum(edge)
