# geneticTSP
A genetic algorithms approach to solve the travel salesman problem

## Case studies
### Djibouti (38 cities)

Source: <http://www.math.uwaterloo.ca/tsp/world/djtour.html>

* Best result using GA: 6659.431532931466 (25 generations, kept until generation 120, then stopped)
* Optimal result: 6656
* Error: 0.05 %

Note: the difference in the first result can be caused by numerical errors (In the TSPLIB norm, the travel cost between each pair of cities is the Euclidean distance between the points rounded to the nearest integer (not the distance rounded to two decimal places). Source: <https://cs.stackexchange.com/questions/56597/tsp-problem-with-a-benchmark-data>). Since the fitness is calculated as 1/distance, this can create some differences.

![Djibouti path using GA](https://raw.githubusercontent.com/nicolopinci/geneticTSP/master/img/djibouti.png "Djibouti path using GA")

<details>
  <summary>Show the results</summary>

```
1 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
2 generations (distance: 6770.076921715761 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13])
3 generations (distance: 6770.07692171576 and best path: [18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13, 19])
4 generations (distance: 6770.07692171576 and best path: [18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13, 19])
5 generations (distance: 6770.07692171576 and best path: [18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 25, 26, 15, 13, 19])
6 generations (distance: 6762.965724246172 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 26, 25, 15, 13])
7 generations (distance: 6762.965724246172 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 26, 25, 15, 13])
8 generations (distance: 6762.965724246172 and best path: [19, 18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 26, 25, 15, 13])
9 generations (distance: 6762.965724246171 and best path: [18, 17, 16, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 20, 23, 26, 25, 15, 13, 19])
10 generations (distance: 6720.065411404567 and best path: [20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 25, 26])
11 generations (distance: 6720.065411404567 and best path: [20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 25, 26])
12 generations (distance: 6720.065411404567 and best path: [20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 25, 26])
13 generations (distance: 6720.065411404567 and best path: [20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 25, 26])
14 generations (distance: 6720.065411404567 and best path: [20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 23, 25, 26])
15 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
16 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
17 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
18 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
19 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
20 generations (distance: 6663.213801224398 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 12, 11, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
21 generations (distance: 6659.906740386759 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 11, 12, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
22 generations (distance: 6659.906740386759 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 11, 12, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
23 generations (distance: 6659.906740386759 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 11, 12, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
24 generations (distance: 6659.906740386759 and best path: [23, 20, 15, 13, 16, 18, 19, 17, 11, 12, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
25 generations (distance: 6659.431532931466 and best path: [23, 20, 15, 13, 16, 17, 18, 19, 11, 12, 9, 8, 7, 6, 5, 3, 4, 2, 1, 10, 14, 21, 29, 30, 32, 35, 37, 38, 33, 34, 36, 31, 27, 28, 24, 22, 25, 26])
```

</details>




### Western Sahara (29 cities)

Source: <http://www.math.uwaterloo.ca/tsp/world/witour.html>

* Best result using GA: 27601.173774493756 (44 generations)
* Optimal result: 27603
* Error: 0 %

Note: the difference in the result can be caused by numerical errors (In the TSPLIB norm, the travel cost between each pair of cities is the Euclidean distance between the points rounded to the nearest integer (not the distance rounded to two decimal places). Source: <https://cs.stackexchange.com/questions/56597/tsp-problem-with-a-benchmark-data>. Since the fitness is calculated as 1/distance, this can create some differences.

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
