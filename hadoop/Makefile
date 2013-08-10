HADOOP = ${HADOOP_HOME}/bin/hadoop

APP = multiplefiles
SRC = src/*.java 
OUT = out

$(APP): $(SRC) 
	mkdir -p $(OUT) 
	javac -classpath `$(HADOOP) classpath` -d $(OUT) $(SRC) 
	jar -cvf $(APP).jar -C $(OUT) .

clean: 
	rm -rf $(OUT) *.jar .
