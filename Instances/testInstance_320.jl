edge = [1 2; 1 12; 2 3; 2 5; 2 8; 2 9; 2 10; 2 13; 2 14; 3 8; 3 10; 3 12; 4 2; 4 3; 4 6; 4 7; 4 8; 4 10; 4 11; 4 13; 4 15; 5 3; 5 7; 5 9; 5 12; 5 13; 5 14; 6 2; 6 10; 6 13; 6 15; 7 2; 7 3; 7 8; 7 9; 7 11; 8 4; 8 6; 8 11; 8 13; 9 2; 9 8; 9 10; 9 14; 10 2; 10 3; 10 11; 10 14; 10 15; 11 4; 11 5; 11 8; 11 12; 11 13; 11 15; 12 2; 12 5; 12 6; 12 7; 12 10; 12 13; 12 14; 13 6; 13 12; 13 14; 13 15; 14 4; 14 5; 14 8; 14 10; 14 15]
cL_orig = [2.0, 470.0, 7.0, 94.0, 191.0, 219.0, 186.0, 68.0, 290.0, 39.0, 59.0, 118.0, 67.0, 25.0, 18.0, 2.0, 140.0, 118.0, 93.0, 26.0, 189.0, 13.0, 0.0, 29.0, 75.0, 172.0, 133.0, 88.0, 102.0, 169.0, 431.0, 87.0, 28.0, 34.0, 70.0, 170.0, 200.0, 40.0, 43.0, 3.0, 222.0, 41.0, 0.0, 47.0, 52.0, 150.0, 49.0, 197.0, 40.0, 25.0, 118.0, 75.0, 2.0, 32.0, 115.0, 444.0, 140.0, 34.0, 202.0, 51.0, 48.0, 1.0, 27.0, 17.0, 22.0, 82.0, 249.0, 436.0, 64.0, 104.0, 38.0]
cU_orig = [42.0, 470.0, 7.0, 102.0, 191.0, 341.0, 186.0, 68.0, 290.0, 39.0, 527.0, 118.0, 67.0, 25.0, 30.0, 2.0, 140.0, 118.0, 93.0, 26.0, 189.0, 25.0, 6.0, 29.0, 75.0, 172.0, 249.0, 88.0, 102.0, 169.0, 431.0, 87.0, 28.0, 34.0, 70.0, 170.0, 200.0, 40.0, 43.0, 3.0, 222.0, 41.0, 4.0, 251.0, 52.0, 276.0, 49.0, 197.0, 40.0, 25.0, 118.0, 77.0, 2.0, 32.0, 115.0, 444.0, 140.0, 252.0, 202.0, 51.0, 48.0, 11.0, 27.0, 17.0, 22.0, 82.0, 249.0, 436.0, 70.0, 104.0, 38.0]
d = [10.0, 35.0, 32.0, 35.0, 36.0, 33.0, 13.0, 2.0, 3.0, 9.0, 11.0, 18.0, 1.0, 8.0, 25.0, 28.0, 27.0, 6.0, 18.0, 18.0, 37.0, 12.0, 47.0, 42.0, 14.0, 40.0, 37.0, 41.0, 27.0, 25.0, 22.0, 44.0, 5.0, 12.0, 21.0, 27.0, 2.0, 7.0, 7.0, 24.0, 31.0, 5.0, 24.0, 28.0, 2.0, 29.0, 47.0, 21.0, 27.0, 47.0, 13.0, 27.0, 23.0, 22.0, 15.0, 35.0, 4.0, 37.0, 2.0, 17.0, 24.0, 39.0, 11.0, 47.0, 34.0, 19.0, 43.0, 38.0, 50.0, 2.0, 30.0]
Len = length(d)

yy = [1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]

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
