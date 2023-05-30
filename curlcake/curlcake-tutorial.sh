cat curlcake-examples
cd curlcake/
curlcake --help

#4mer
java -jar curlcake.jar 4 28 1 100 incomplete.txt complete.txt ./RNAshapes  0
java -jar curlcake.jar 4 24 1 100 incomplete.txt complete24s0.txt ./RNAshapes  0
java -jar curlcake.jar 4 40 1 100 incomplete.txt complete4-40s0.txt ./RNAshapes  0
java -jar curlcake.jar 4 20 1 100 incomplete.txt complete4-20s0.txt ./RNAshapes  0
java -jar curlcake.jar 4 40 1 100 incomplete-4-40-1-100.txt complete-4-40-1-100-s1028.txt ./RNAshapes 1028  

cd curlcake/
java -jar curlcake.jar 3 48 1 100 incomplete.txt complete-3-48-s1028.txt ./RNAshapes  1028
java -jar curlcake.jar 3 48 1 100 incomplete.txt complete-3-48-s102.txt ./RNAshapes  102
java -jar curlcake.jar 3 51 1 100 incomplete.txt complete-3-51-s100.txt ./RNAshapes  100
java -jar curlcake.jar 3 48 1 100 incomplete.txt complete-3-48-s777.txt ./RNAshapes 777
java -jar curlcake.jar 3 30 1 100 incomplete.txt complete-3-48-s777.txt ./RNAshapes 1028
java -jar curlcake.jar 5 40 1 100 5.40.1.100.s1028_incomplete.txt 5.40.1.100.s1028_complete.txt ./RNAshapes 1028
java -jar curlcake.jar 5 40 1 100 5.40.1.100.s456_incomplete.txt 5.40.1.100.s456_complete.txt ./RNAshapes 456

#content of curlcake example file
java -jar curlcake.jar 5 50 1 100 incomplete_5_50_1_100.txt complete_5_50_1_100.txt ./RNAshapes 0
java -jar curlcake.jar 5 25 1 100 incomplete_5_25_1_100.txt complete_5_25_1_100.txt ./RNAshapes 0
java -jar curlcake.jar 5 30 1 100 incomplete_5_30_1_100.txt complete_5_30_1_100.txt ./RNAshapes 0
java -jar curlcake.jar 5 40 1 100 incomplete_5_40_1_100.txt complete_5_40_1_100.txt ./RNAshapes 0

java -jar curlcake.jar 5 35 1 100 incomplete_5_35.txt complete_5_35.txt ./RNAshapes 0

java -jar curlcake.jar 6 40 1 100 incompleteS456.txt complete.6.40.1.100.s456.txt ./RNAshapes 456 
java -jar curlcake.jar 7 40 1 100 incomplete7.S456.txt complete.7.40.1.100.S456.txt ./RNAshapes 456 

 
