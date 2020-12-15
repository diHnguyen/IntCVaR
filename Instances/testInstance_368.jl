edge = [1 4; 1 8; 1 14; 1 15; 2 3; 2 5; 2 10; 2 13; 3 14; 4 3; 4 6; 4 8; 4 9; 4 11; 4 14; 4 15; 5 2; 5 3; 5 6; 5 8; 5 9; 5 11; 6 4; 6 9; 6 12; 6 13; 6 14; 6 15; 7 8; 7 10; 7 12; 7 13; 7 14; 7 15; 8 2; 8 5; 8 7; 8 13; 8 14; 9 2; 9 4; 9 6; 9 10; 9 12; 9 14; 10 4; 10 5; 10 6; 10 8; 10 11; 10 12; 11 4; 11 12; 11 13; 11 14; 12 2; 12 5; 12 8; 12 10; 12 13; 12 14; 13 5; 13 7; 13 11; 13 12; 13 14; 13 15; 14 4; 14 11]
cL_orig = [134.0, 195.0, 427.0, 264.0, 21.0, 118.0, 399.0, 244.0, 5.0, 21.0, 60.0, 37.0, 213.0, 159.0, 61.0, 376.0, 71.0, 35.0, 47.0, 6.0, 56.0, 243.0, 34.0, 63.0, 124.0, 137.0, 38.0, 284.0, 19.0, 118.0, 40.0, 33.0, 222.0, 33.0, 158.0, 61.0, 25.0, 92.0, 136.0, 227.0, 148.0, 17.0, 26.0, 46.0, 15.0, 131.0, 199.0, 22.0, 12.0, 2.0, 55.0, 132.0, 40.0, 6.0, 82.0, 308.0, 258.0, 51.0, 1.0, 9.0, 26.0, 284.0, 182.0, 15.0, 30.0, 45.0, 69.0, 27.0, 51.0]
cU_orig = [134.0, 279.0, 427.0, 264.0, 21.0, 118.0, 399.0, 244.0, 209.0, 21.0, 60.0, 37.0, 213.0, 159.0, 61.0, 376.0, 71.0, 35.0, 47.0, 12.0, 56.0, 243.0, 34.0, 115.0, 392.0, 137.0, 38.0, 284.0, 19.0, 118.0, 398.0, 33.0, 222.0, 101.0, 158.0, 61.0, 25.0, 92.0, 136.0, 227.0, 148.0, 17.0, 26.0, 46.0, 15.0, 131.0, 199.0, 22.0, 12.0, 2.0, 55.0, 132.0, 40.0, 6.0, 82.0, 308.0, 258.0, 339.0, 1.0, 69.0, 26.0, 284.0, 182.0, 181.0, 30.0, 45.0, 69.0, 53.0, 97.0]
d = [5.0, 45.0, 22.0, 4.0, 40.0, 49.0, 3.0, 17.0, 15.0, 38.0, 40.0, 14.0, 27.0, 19.0, 23.0, 13.0, 32.0, 3.0, 44.0, 32.0, 25.0, 18.0, 25.0, 16.0, 42.0, 27.0, 14.0, 3.0, 49.0, 43.0, 36.0, 2.0, 48.0, 10.0, 14.0, 9.0, 1.0, 38.0, 9.0, 36.0, 28.0, 20.0, 24.0, 42.0, 35.0, 50.0, 12.0, 28.0, 10.0, 25.0, 22.0, 19.0, 33.0, 47.0, 22.0, 48.0, 40.0, 40.0, 22.0, 20.0, 45.0, 2.0, 40.0, 20.0, 36.0, 44.0, 28.0, 14.0, 45.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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