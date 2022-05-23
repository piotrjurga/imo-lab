main: main.cpp
	g++ -g main.cpp -o main -std=c++17

release: main.cpp
	g++ -O3 main.cpp -o main -std=c++17
