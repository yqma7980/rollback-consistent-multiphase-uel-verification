*Node
1, 0.0, 0.0
2, 1.0, 0.0
3, 2.0, 0.0
4, 0.0, 1.0
5, 1.0, 1.0
6, 2.0, 1.0
7, 0.0, 2.0
8, 1.0, 2.0
9, 2.0, 2.0
*User Element, Type=U1, Nodes=4, Coordinates=2, Properties=27, I Properties=1, Variables=192, Unsymm
1, 2, 8, 14
*Element, Type=U1, Elset=PATCH
10287, 1, 2, 5, 4
10288, 2, 3, 6, 5
10289, 4, 5, 8, 7
10290, 5, 6, 9, 8
*UEL Property, Elset=PATCH
1.0E9, 0.25, 0.80, 1.0E-9, 1.0E-3, 700.0, 1.0E-3, 0.20
0.05, 1.0E4, 2.0, 2.0, 2.0, 5.0E-3, 0.10, 1.0E-14
0.0, 0.0, 0.0, 0.0, 1.0E-20, 1.0E-8, 0.0, 0.0
0.0, 0.0, 0.0
0
*Nset, Nset=ALLNODES, Generate
1, 9, 1
*Nset, Nset=LEFT
1, 4, 7
*Nset, Nset=RIGHT
3, 6, 9
*Nset, Nset=BOTTOM
1, 2, 3
