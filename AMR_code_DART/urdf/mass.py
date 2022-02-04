fopen = open("anymal.urdf","r",encoding="utf-8")
mass = 0
cont=0
for i in fopen:
    if "mass" in i:
        cont+=1
        l = i.split(" ")
        for j in l:
            if "value" in j:
                l1 = j.split("value=")
                l2 = l1[1][1:-3]
                mass+=float(l2)
fopen.close()
print(mass,cont)
