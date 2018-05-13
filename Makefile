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
	echo 'data="$$(java -jar tester.jar -exec ./log/${timestamp}.bin -debug -seed $$1 | tee /dev/stderr | grep '\''{"seed":'\'')" ; flock log/${timestamp}.lock echo "$$data" >> log/${timestamp}.json' > log/${timestamp}.sh
	parallel --jobs 50% -- bash log/${timestamp}.sh {} ::: $$(seq 100)
	cat log/${timestamp}.json | python3 stat-results.py table
	cat log/${timestamp}.json | python3 stat-results.py summary
