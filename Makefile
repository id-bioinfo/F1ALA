F1ALA: TIPars.java
	javac -classpath .:lib/beast.jar:lib/gson-2.8.5.jar TIPars.java;\
	jar cvfm F1ALA.jar MANIFEST.MF *class;\
	rm -rf *class
