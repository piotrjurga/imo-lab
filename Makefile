main: main.cpp
	g++ -g main.cpp -o main

release: main.cpp
	g++ -O3 main.cpp -o main
