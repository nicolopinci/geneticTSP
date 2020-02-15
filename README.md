# geneticTSP
A genetic algorithms approach to solve the travel salesman problem

## Case studies
### Djibouti (38 cities)

Source: <http://www.math.uwaterloo.ca/tsp/world/djtour.html>

* Best result using GA: 6663.213801224396 (47 generations)
* Optimal result: 6656
* Error: 0.108 %

Note: this result refers to the case in which the greedy chromosomes are create at the beginning of the evolutionary process.

![Djibouti path using GA](https://raw.githubusercontent.com/nicolopinci/geneticTSP/master/img/djibouti.png "Djibouti path using GA")

<details>
  <summary>Show the results</summary>

```
1 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
2 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
3 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
4 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
5 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
6 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
7 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
8 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
9 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
10 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
11 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
12 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
13 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
14 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
15 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
16 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
17 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
18 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
19 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
20 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
21 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
22 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
23 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
24 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
25 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
26 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
27 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
28 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
29 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
30 generations (distance: 6769.69757480714 and best path: [16, 19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
31 generations (distance: 6715.63024386395 and best path: [8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 19, 18, 17, 16, 12, 11, 9])
32 generations (distance: 6715.63024386395 and best path: [8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 19, 18, 17, 16, 12, 11, 9])
33 generations (distance: 6715.63024386395 and best path: [8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 19, 18, 17, 16, 12, 11, 9])
34 generations (distance: 6715.63024386395 and best path: [8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 19, 18, 17, 16, 12, 11, 9])
35 generations (distance: 6715.63024386395 and best path: [8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 19, 18, 17, 16, 12, 11, 9])
36 generations (distance: 6715.250896955328 and best path: [19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
37 generations (distance: 6715.250896955328 and best path: [19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
38 generations (distance: 6715.250896955328 and best path: [19, 18, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
39 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
40 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
41 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
42 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
43 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
44 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
45 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
46 generations (distance: 6714.730474523721 and best path: [18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 26, 25, 20, 15, 13, 16])
47 generations (distance: 6663.213801224396 and best path: [19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26, 23, 20, 15, 13, 16, 18])
```

</details>




### Western Sahara (29 cities)

Source: <http://www.math.uwaterloo.ca/tsp/world/djtour.html>

* Best result using GA: 27601.173774493756 (44 generations)
* Optimal result: 27603
* Error: 0 %

Note: the difference in the result can be caused by numerical errors

![Sahara path using GA](https://raw.githubusercontent.com/nicolopinci/geneticTSP/master/img/sahara.png "Sahara path using GA")

<details>
  <summary>Show the results</summary>

```
1 generations (distance: 32161.402974904493 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 14, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 15, 9, 7, 3])
2 generations (distance: 32161.402974904493 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 14, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 15, 9, 7, 3])
3 generations (distance: 31851.693014269145 and best path: [5, 9, 7, 3, 4, 8, 12, 13, 14, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 15, 11, 10, 6, 2, 1])
4 generations (distance: 31851.693014269145 and best path: [5, 9, 7, 3, 4, 8, 12, 13, 14, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 15, 11, 10, 6, 2, 1])
5 generations (distance: 31851.693014269145 and best path: [5, 9, 7, 3, 4, 8, 12, 13, 14, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 15, 11, 10, 6, 2, 1])
6 generations (distance: 30772.021521374416 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 14, 9, 7, 3])
7 generations (distance: 30770.43704665055 and best path: [8, 4, 5, 6, 2, 1, 11, 10, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 25, 27, 24, 14, 9, 7, 3])
8 generations (distance: 30521.0077559403 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7, 3])
9 generations (distance: 30521.0077559403 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7, 3])
10 generations (distance: 30521.0077559403 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7, 3])
11 generations (distance: 30521.0077559403 and best path: [8, 4, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7, 3])
12 generations (distance: 30249.81957574448 and best path: [3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7])
13 generations (distance: 29979.343574930048 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
14 generations (distance: 29979.343574930048 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
15 generations (distance: 29979.343574930048 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
16 generations (distance: 29624.628506282268 and best path: [3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7])
17 generations (distance: 29624.628506282268 and best path: [3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 16, 24, 27, 25, 14, 9, 7])
18 generations (distance: 29548.4573210376 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
19 generations (distance: 29354.152505467835 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
20 generations (distance: 29354.152505467835 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
21 generations (distance: 29354.152505467835 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
22 generations (distance: 29354.152505467835 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 24, 27, 16, 14, 9])
23 generations (distance: 28923.266251575387 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
24 generations (distance: 28923.266251575387 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
25 generations (distance: 28923.266251575387 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
26 generations (distance: 28923.266251575387 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
27 generations (distance: 28923.266251575387 and best path: [7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 13, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14, 9])
28 generations (distance: 28364.57501194755 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
29 generations (distance: 28364.57501194755 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 17, 18, 19, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
30 generations (distance: 27936.376114029186 and best path: [14, 13, 9, 7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16])
31 generations (distance: 27936.376114029186 and best path: [14, 13, 9, 7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16])
32 generations (distance: 27936.376114029186 and best path: [14, 13, 9, 7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16])
33 generations (distance: 27936.376114029186 and best path: [14, 13, 9, 7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16])
34 generations (distance: 27936.376114029186 and best path: [14, 13, 9, 7, 3, 4, 8, 5, 6, 2, 1, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16])
35 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
36 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
37 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
38 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
39 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
40 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
41 generations (distance: 27739.383942485336 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 22, 23, 21, 29, 28, 26, 20, 25, 27, 24, 16, 14])
42 generations (distance: 27614.4780817462 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 21, 23, 22, 29, 28, 26, 20, 25, 27, 24, 16, 14])
43 generations (distance: 27614.4780817462 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 21, 23, 22, 29, 28, 26, 20, 25, 27, 24, 16, 14])
44 generations (distance: 27601.173774493756 and best path: [13, 9, 7, 3, 4, 8, 5, 1, 2, 6, 10, 11, 12, 15, 19, 18, 17, 21, 22, 23, 29, 28, 26, 20, 25, 27, 24, 16, 14])
```
</details>
