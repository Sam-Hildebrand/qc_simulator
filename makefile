all:
	g++ main.cpp -o qc_simulator

clean:
	rm -f qc_simulator

run: all
	./qc_simulator