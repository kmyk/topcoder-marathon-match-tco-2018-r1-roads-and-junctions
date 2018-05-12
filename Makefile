.PHONY: build
.DEFAULT: build

PROBLEM := RoadsAndJunctions
CXX := g++
CXXFLAGS := -std=c++11 -Wall -O2 -g -DLOCAL

build: a.out tester.jar

run: a.out tester.jar
	java -jar tester.jar -exec ./a.out

a.out: main.cpp ${PROBLEM}.cpp
	${CXX} ${CXXFLAGS} $<

tester.jar: ${PROBLEM}Vis.java
	javac $<
	jar cvfe tester.jar ${PROBLEM}Vis *.class

URL := https://community.topcoder.com/longcontest/?module=ViewProblemStatement&compid=63555&rd=17153
submit:
	oj submit '${URL}' --language C++ ${PROBLEM}.cpp -y --open
submit/full:
	oj submit '${URL}' --language C++ ${PROBLEM}.cpp -y --open --full

standings:
	oj get-standings '${URL}' --format=csv | \
		sed 's/\w\w\.\w\w\.\w\w\w\w \w\w:\w\w:\w\w/& EST/ ; y/ /~/ ; :1 ; s/,,/,-,/ ; t1 ; s/$$/,|/' | \
		column -t -s , | \
		sed 's/\( \+\)|\?/\1| /g ; s/^/| / ; s/.$$// ; y/~/ / ; 1 { p ; s/[^|]/-/g }'

timestamp := $(shell date +%s)
score: a.out tester.jar
	-mkdir log
	cp a.out log/${timestamp}.bin
	for seed in $$(seq 1 100) ; do \
		java -jar tester.jar -exec ./log/${timestamp}.bin -debug -seed $$seed | \
			tee /dev/stderr | \
			grep '{"seed":' >> log/${timestamp}.json ; \
	done
	cat log/${timestamp}.json | python3 stat-results.py table
	cat log/${timestamp}.json | python3 stat-results.py summary

plot-gradient:
	${CXX} ${CXXFLAGS} plot-gradient.cpp -o plot-gradient.bin
	for seed in $$(seq 1 100) ; do \
		file test/$$seed.in ; \
		N=1000 time ./plot-gradient.bin < test/$$seed.in | python3 plot-gradient.py test/$$seed.in /dev/stdin --save $$seed.png --no-frame ; \
		file $$seed.png ; \
	done
