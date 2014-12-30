# for linux
cd ..
rm -f kill.txt
rm -f *.class
javac -Xlint:all -Xmaxerrs 5 -cp .:lib/* Analysis.java


if [ $? -eq 0 ]; then
    echo Compiled.
    
    java -Xmx30g -XX:ParallelGCThreads=8 -cp .:lib/* Analysis $1
fi

rm -f *.class
