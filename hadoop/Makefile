SHELL = /bin/bash
APP = multiplefiles
HADOOP=/Users/eterna/hadoop-2.4.0
SRC = src/*.java
EBIRD = relevantelephant/*.java

$(APP): $(SRC) 
	rm -rf out
	mkdir -p out
	mkdir -p  el
	javac -classpath $(shell find $(HADOOP) -name "*.jar" | tr '\n' ':' ) -d out $(SRC)
	cp $(HADOOP)/share/hadoop/tools/lib/hadoop-streaming*.jar out/
	cd out; mv hadoop-streaming*.jar hadoop-streaming-mod.jar; jar uf hadoop-streaming*.jar org/apache/hadoop/streaming/* edu/jhu/cs/MultipleOutputs*
	
clean: 
	rm -rf out *.jar .
