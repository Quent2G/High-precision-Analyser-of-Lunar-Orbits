

with open("Utils/EGM96.gfc") as text:
    Lines = text.readlines()

with open("Utils/C_EGM96.gfc","w") as NT:
    l = 1
    for line in Lines:
        if l==1:
            NT.write("0   0   1.0   0.0\n")
        else:
            data = line.split()
            NT.write(f"{data[1]}   {data[2]}   {data[3]}   {data[4]}\n")
        l+=1