main: main.cpp
	mkdir -p bin &
	g++ main.cpp -lgmp -Ofast -o bin/main

run: main.cpp
	mkdir -p bin &
	g++ main.cpp -lgmp -Ofast -o bin/main && bin/main $(d)
