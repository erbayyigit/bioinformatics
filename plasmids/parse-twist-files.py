import csv 
'''
  input file is Twist order file
'''
twist_file = open('TadA8.20_clones.csv')
csv_reader = csv.DictReader(twist_file)
for row in csv_reader:
    mystring =f">{row['Name']}\n{row['Construct sequence']}"
    print(mystring)
    
