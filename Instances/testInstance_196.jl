edge = [1 7; 1 14; 2 5; 2 6; 2 7; 2 12; 2 13; 2 15; 3 4; 3 7; 3 14; 3 15; 4 2; 4 11; 4 12; 4 13; 4 14; 5 10; 5 12; 5 14; 5 15; 6 4; 6 11; 6 12; 7 9; 7 10; 7 11; 7 12; 7 13; 8 9; 8 13; 9 2; 9 4; 9 13; 10 2; 10 4; 10 6; 10 7; 10 12; 10 15; 11 4; 11 5; 11 14; 12 3; 12 5; 12 10; 12 15; 13 3; 13 4; 13 6; 13 7; 13 11; 13 15; 14 3; 14 5; 14 10]
cL_orig = [258.0, 81.0, 77.0, 138.0, 162.0, 362.0, 200.0, 120.0, 40.0, 22.0, 144.0, 252.0, 36.0, 322.0, 133.0, 199.0, 293.0, 157.0, 188.0, 72.0, 265.0, 8.0, 141.0, 119.0, 59.0, 129.0, 192.0, 218.0, 185.0, 23.0, 79.0, 165.0, 110.0, 58.0, 89.0, 109.0, 158.0, 29.0, 75.0, 135.0, 197.0, 200.0, 127.0, 17.0, 78.0, 67.0, 63.0, 471.0, 343.0, 20.0, 224.0, 10.0, 41.0, 354.0, 11.0, 9.0]
cU_orig = [258.0, 587.0, 77.0, 138.0, 162.0, 512.0, 200.0, 696.0, 40.0, 138.0, 144.0, 252.0, 162.0, 322.0, 133.0, 199.0, 293.0, 157.0, 188.0, 72.0, 265.0, 8.0, 141.0, 119.0, 59.0, 137.0, 192.0, 218.0, 185.0, 23.0, 79.0, 165.0, 184.0, 58.0, 89.0, 109.0, 158.0, 91.0, 75.0, 135.0, 197.0, 200.0, 127.0, 17.0, 400.0, 67.0, 217.0, 471.0, 343.0, 332.0, 224.0, 10.0, 41.0, 354.0, 11.0, 9.0]
d = [25.0, 32.0, 7.0, 20.0, 47.0, 1.0, 2.0, 44.0, 28.0, 31.0, 1.0, 37.0, 45.0, 8.0, 41.0, 44.0, 18.0, 23.0, 30.0, 39.0, 18.0, 46.0, 49.0, 24.0, 37.0, 21.0, 30.0, 21.0, 16.0, 40.0, 9.0, 27.0, 34.0, 42.0, 34.0, 13.0, 12.0, 24.0, 27.0, 28.0, 1.0, 1.0, 26.0, 47.0, 23.0, 50.0, 41.0, 27.0, 31.0, 19.0, 9.0, 44.0, 7.0, 39.0, 16.0, 48.0]
Len = length(d)

yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]

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
