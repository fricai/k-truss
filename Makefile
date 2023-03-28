main:
	mpic++ -Wall -std=c++17 -O3 -funroll-loops code/main.cpp -o a3

clean:
	rm a3
