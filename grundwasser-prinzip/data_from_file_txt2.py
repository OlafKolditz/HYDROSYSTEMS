import matplotlib.pyplot as plt
#import csv

#with open('out.txt') as f:
#    contents = f.read()
#    print(contents)
	
#f.close()

fig, ax = plt.subplots(ncols=2, figsize=(12,5))

x = []
h = []
i = 0
with open("out1.txt", "r") as text_file:
    for line in text_file:
        xp, hp = line.strip().split(", ")
        i=i+1
        x.append(xp)
        h.append(float(hp))
    ax[0].plot(x,h, label='i='+str(i/3))
		
print(h)

ax[0].set_title('Hydrosystemanalyse\nExercise Prinzipbeispiel Grundwasser')
ax[0].set_xlabel('x')
ax[0].set_ylabel('h')
ax[0].legend()
ax[0].grid()

x = []
h = []
i = 0
with open("out2.txt", "r") as text_file:
    for line in text_file:
        xp, hp = line.strip().split(", ")
        i=i+1
        x.append(xp)
        h.append(float(hp))
    ax[1].plot(x,h, label='i='+str(i/3))

print(h)

ax[1].set_title('Hydrosystemanalyse\nExercise Prinzipbeispiel Grundwasser')
ax[1].set_xlabel('x')
ax[1].set_ylabel('h')
ax[1].legend()
ax[1].grid()

plt.title('Hydrosystemanalyse\nExercise Prinzipbeispiel Grundwasser')
plt.savefig("gauss-seidel-plt.png")
plt.show()
