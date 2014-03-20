APP = multiplefiles
HADOOP=/bin/hadoop
SRC = src/*.java 
OUT = out

$(APP): $(SRC) 
	mkdir -p $(OUT) 
	javac -classpath $(HADOOP)/hadoop-core-1.0.3.jar:$(HADOOP)/lib/commons-cli-1.2.jar -d $(OUT) $(SRC)
	jar -cvf $(APP).jar -C $(OUT) .

clean: 
	rm -rf $(OUT) *.jar .
