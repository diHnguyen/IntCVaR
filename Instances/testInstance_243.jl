edge = [1 3; 1 5; 1 7; 1 12; 1 13; 2 3; 2 6; 2 7; 2 14; 3 2; 3 4; 3 5; 3 6; 3 8; 4 2; 4 6; 4 8; 4 11; 4 12; 4 13; 4 14; 4 15; 5 6; 5 7; 5 9; 5 11; 5 13; 6 2; 6 3; 6 7; 6 10; 6 14; 7 2; 7 4; 7 5; 7 8; 7 10; 7 15; 8 10; 8 15; 9 6; 9 12; 9 13; 9 14; 10 2; 10 3; 10 8; 10 11; 10 14; 11 2; 11 12; 11 13; 11 14; 12 3; 12 5; 12 8; 12 9; 12 11; 13 3; 13 5; 14 2; 14 5; 14 6; 14 7; 14 8; 14 11; 14 12]
cL_orig = [10.0, 89.0, 132.0, 424.0, 20.0, 22.0, 39.0, 23.0, 214.0, 15.0, 3.0, 15.0, 61.0, 171.0, 77.0, 8.0, 55.0, 51.0, 176.0, 222.0, 65.0, 76.0, 10.0, 54.0, 177.0, 297.0, 192.0, 135.0, 57.0, 17.0, 32.0, 202.0, 195.0, 144.0, 50.0, 35.0, 10.0, 318.0, 85.0, 192.0, 9.0, 31.0, 135.0, 21.0, 301.0, 243.0, 52.0, 23.0, 152.0, 159.0, 5.0, 25.0, 123.0, 197.0, 214.0, 46.0, 2.0, 32.0, 194.0, 145.0, 507.0, 101.0, 130.0, 202.0, 52.0, 142.0, 56.0]
cU_orig = [150.0, 89.0, 132.0, 424.0, 20.0, 32.0, 253.0, 23.0, 214.0, 15.0, 3.0, 15.0, 61.0, 171.0, 77.0, 32.0, 63.0, 51.0, 176.0, 222.0, 191.0, 76.0, 10.0, 54.0, 177.0, 297.0, 282.0, 135.0, 57.0, 17.0, 32.0, 202.0, 195.0, 144.0, 50.0, 35.0, 10.0, 318.0, 85.0, 192.0, 73.0, 151.0, 135.0, 361.0, 301.0, 243.0, 52.0, 23.0, 152.0, 159.0, 5.0, 25.0, 123.0, 537.0, 214.0, 46.0, 2.0, 32.0, 194.0, 145.0, 507.0, 101.0, 262.0, 202.0, 52.0, 142.0, 56.0]
d = [27.0, 32.0, 20.0, 40.0, 26.0, 48.0, 33.0, 28.0, 4.0, 13.0, 31.0, 21.0, 17.0, 18.0, 29.0, 16.0, 6.0, 47.0, 37.0, 36.0, 43.0, 38.0, 3.0, 40.0, 38.0, 36.0, 18.0, 41.0, 12.0, 23.0, 26.0, 9.0, 2.0, 21.0, 45.0, 4.0, 42.0, 19.0, 37.0, 29.0, 48.0, 48.0, 48.0, 27.0, 5.0, 28.0, 40.0, 2.0, 2.0, 22.0, 36.0, 23.0, 24.0, 18.0, 6.0, 28.0, 33.0, 22.0, 10.0, 47.0, 44.0, 26.0, 6.0, 19.0, 3.0, 19.0, 27.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

c_orig = 0.5*(cL_orig+cU_orig)

SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)


p = [1.0]

g = [SP_init]

h = [0.0]


origin = 1

destination =15

last_node = maximum(edge)
all_nodes = collect(1:last_node)

M_orig = zeros(Len)

for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end

case = 0
delta1 = 1e-6
delta2 = 5
last_node = maximum(edge)