

with open("input/model2convert/Earth_EGM2008.txt") as text:
    Lines = text.readlines()

with open("input/model2convert/Earth_C-EGM2008.txt","w") as NT:
    l = 1
    for line in Lines:
        if l==1:
            NT.write("0   0   1.0   0.0\n")
        else:
            data = line.split()
            NT.write(f"{data[0]}   {data[1]}   {data[2]}   {data[3]}\n")
        l+=1