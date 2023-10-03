TIPars2: TIPars.java
	javac -classpath .:lib/beast.jar:lib/gson-2.8.5.jar TIPars.java;\
	jar cvfm TIPars2.jar MANIFEST.MF *class;\
	rm -rf *class
