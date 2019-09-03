def writeStatus(s):
    if s == "c":
        with open('cNMRFile.txt', 'w') as text:
            text.write("C")
    elif s == "f":
        with open('cNMRFile.txt', 'w') as text:
            text.write("F")
    print('Status Updated Succesfully')

#Determine type of NMR data - txt communication with MatLab
with open('cNMRFile.txt') as f:
    type = f.readline()

#Topspin Processing
if type == "TOPSPIN":
    print("TopSpin (Bruker) Data Detected...")
    inputFile = "NMR_TOPSPIN.txt"
    file = open(inputFile, 'r')
    bigList = [[] for i in range(0,2048)]
    counter = 0
    column = 0

    for i in file:
        if i[0] != "#":
            i = i.strip('\n')
            bigList[counter].append(i)
            counter += 1
        else:
            if counter == 2048:
                counter = 0
    
    fileOut = open('nmrData.dat', 'w')
    for i in range(0, len(bigList)):
        for j in range(0, len(bigList[i])):
            fileOut.write(bigList[i][j])
            fileOut.write(",")
        fileOut.write('\n')
    fileOut.close()
    writeStatus("c")
    print("TopSpin NMR data converted to MATLAB compatible delimited data.")

#Agilent Processing
elif type == "AGILENT":
    writeStatus("f")
    print("Agilent Data Detected...")
    print("This type of data is not yet supported!")

#If data is already processed
elif type == "DELIM":
    writeStatus("c")
    print("Delimited Data Selected, Python Not Required")

#If data identity can't be determined, or not known type
else:
    writeStatus("f")
    print("Either this wasn't launched by MATLAB, or the data is crap!")